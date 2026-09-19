import pandas as pd
import pytest

from eacbp.artifact.registry import ArtifactRegistry
from eacbp.auditor import ScientificAuditor
from eacbp.knowledge.capability import KnowledgeRetrievalCapability
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus


@pytest.mark.parametrize('prior', [False, True])
def test_real_knowledge_output_has_required_audit_in_both_modes(tmp_path, prior):
    registry = ArtifactRegistry(str(tmp_path))
    method = 'knowledge_engine_prior_v1' if prior else 'knowledge_engine_discovery_v1'
    task = TaskContract(task_id='knowledge', capability='knowledge_retrieval', method=method,
                        parameters={'study_id': 's', 'species': 'human', 'tissue': 'brain',
                                    'prior_guided': prior, 'hypotheses': ['Apoe pathway'], 'target_genes': ['Apoe']},
                        expected_outputs=['table://s/knowledge/v1', 'json://s/report/v1'],
                        validation_requirements=['epistemic_tagging_check'])
    result = KnowledgeRetrievalCapability(method).execute(task, registry)
    report = ScientificAuditor().audit_task(task, result, registry)
    assert report.overall_passed, report.model_dump()
    assert any(check.check_name == 'epistemic_tagging_check' and check.passed for check in report.checks)


@pytest.mark.parametrize('prior', [False, True])
def test_metrics_cannot_substitute_for_persisted_knowledge_labels(tmp_path, prior):
    registry = ArtifactRegistry(str(tmp_path))
    registry.register('table://s/knowledge/v1', pd.DataFrame({'category': ['None'], 'name': ['No evidence']}),
                      ArtifactType.TABLE, 's', 'knowledge', 'fixture')
    registry.register('json://s/report/v1', {'mode': 'prior_guided' if prior else 'discovery',
                      'prior_guided': prior, 'summary': '', 'epistemic_tags': [], 'evidence_nodes': []},
                      ArtifactType.JSON, 's', 'knowledge', 'fixture')
    task = TaskContract(task_id='knowledge', capability='knowledge_retrieval', parameters={'prior_guided': prior},
                        validation_requirements=['epistemic_tagging_check'])
    result = TaskResult(task_id='knowledge', capability='knowledge_retrieval', method_used='fixture',
                        status=TaskStatus.SUCCESS, output_artifacts=['table://s/knowledge/v1', 'json://s/report/v1'],
                        metrics={'epistemic_tags': ['[PRIOR-GUIDED HYPOTHESIS TESTING]', 'mode:unbiased_discovery']})
    report = ScientificAuditor().audit_task(task, result, registry)
    assert not report.overall_passed
    assert any(check.check_name == 'epistemic_tagging_check' and not check.passed for check in report.checks)
