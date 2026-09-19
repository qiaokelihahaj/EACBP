from unittest.mock import patch

import pytest
from pydantic import BaseModel, ConfigDict, Field, ValidationError

from eacbp.artifact.registry import ArtifactRegistry
from eacbp.auditor.base import ValidationReport, ValidationSeverity
from eacbp.capabilities.base import BaseCapability, CapabilityDescriptor
from eacbp.capabilities.registry import CapabilityRegistry
from eacbp.orchestrator.dag import ComputationalDAGPlanner
from eacbp.orchestrator.loop import ScientificOrchestrator
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.evidence import EvidenceNode, EvidenceType
from eacbp.schemas.study import StudyManifest, BiologicalDesign
from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus


class Parameters(BaseModel):
    model_config = ConfigDict(extra='allow')
    count: int = Field(ge=1)


def independent_audit(contract, result, registry):
    payload = registry.get(result.output_artifacts[0])[1]
    report = ValidationReport(auditor_name='independent_fixture', target_task_id=contract.task_id)
    report.add_check('fixture_positive', payload['value'] > 0, ValidationSeverity.ERROR, 'Read persisted value')
    return report


def evidence(contract, result, report, registry):
    return [EvidenceNode(evidence_id=f'E_{contract.task_id}', type=EvidenceType.QC_METRICS,
                         summary='Audited fixture', source_task_id=contract.task_id,
                         source_artifact_uris=result.output_artifacts, audit_passed=True)]


def plan_factory(*, manifest, parameters, tasks):
    return TaskContract(task_id='extension', capability='fixture', parameters=parameters,
                        expected_outputs=[f'json://{manifest.study_id}/fixture/v1'])


class FixtureCapability(BaseCapability):
    def __init__(self, value=1):
        super().__init__('fixture', 'fixture_v1', descriptor=CapabilityDescriptor(
            'fixture', 'fixture_v1', parameter_model=Parameters,
            output_types=[ArtifactType.JSON], required_audit_ids=['fixture_positive'],
            validator=independent_audit, evidence_extractor=evidence, plan_factory=plan_factory))
        self.value = value
        self.calls = 0

    def execute(self, contract, registry):
        self.calls += 1
        assert isinstance(contract.parameters['count'], int)
        uri = contract.expected_outputs[0]
        registry.register(uri, {'value': self.value}, ArtifactType.JSON, 's', contract.task_id, 'fixture')
        return TaskResult(task_id=contract.task_id, capability=self.capability_name,
                          method_used=self.implementation_id, status=TaskStatus.SUCCESS,
                          output_artifacts=[uri])


@pytest.mark.parametrize('value,accepted', [(1, True), (-1, False)])
def test_descriptor_connects_plan_parameters_audit_and_evidence(tmp_path, value, accepted):
    registry = CapabilityRegistry()
    capability = registry.register(FixtureCapability(value))
    orchestrator = ScientificOrchestrator(ArtifactRegistry(str(tmp_path)), registry)
    manifest = StudyManifest(study_id='s', biological_design=BiologicalDesign(species='human', tissue='test'))
    with patch.object(ComputationalDAGPlanner, 'build_study_plan', return_value=[]):
        result = orchestrator.run_study(manifest, {'analysis_extensions': {'fixture': {'count': '2'}}})
    assert capability.calls == 1, result
    assert (result['status'] == 'success') is accepted, result
    assert bool(orchestrator.evidence_graph.evidence_nodes) is accepted
    assert any(check.check_name == 'fixture_positive' for report in orchestrator.audit_reports for check in report.checks)


def test_descriptor_rejects_bad_parameters_and_output_types_before_execution(tmp_path):
    registry = CapabilityRegistry()
    capability = registry.register(FixtureCapability())
    contract = TaskContract(task_id='t', capability='fixture', parameters={'count': 0})
    with pytest.raises(ValidationError):
        registry.prepare_contract(contract)
    contract.parameters = {'count': '2', 'target_branch': 'a'}
    contract.expected_outputs = ['table://s/wrong/v1']
    with pytest.raises(ValueError, match='output types'):
        registry.prepare_contract(contract)
    assert capability.calls == 0
    contract.expected_outputs = ['json://s/right/v1']
    registry.prepare_contract(contract)
    assert contract.parameters == {'count': 2, 'target_branch': 'a'}
    assert contract.validation_requirements == ['fixture_positive']


def test_descriptor_rejects_cross_task_evidence(tmp_path):
    registry = CapabilityRegistry()
    capability = registry.register(FixtureCapability())
    artifacts = ArtifactRegistry(str(tmp_path))
    contract = TaskContract(task_id='t', capability='fixture', parameters={'count': 1},
                            expected_outputs=['json://s/right/v1'])
    result = capability.execute(contract, artifacts)
    report = independent_audit(contract, result, artifacts)
    def wrong_source(*args):
        node = evidence(*args)[0]
        node.source_task_id = 'other'
        return [node]
    capability.descriptor.evidence_extractor = wrong_source
    with pytest.raises(ValueError):
        registry.extract_evidence(contract, result, report, artifacts)
