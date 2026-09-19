from eacbp.auditor.requirements import missing_required_checks
from eacbp.auditor.base import ValidationCheck


def check(name):
    return ValidationCheck(check_name=name, passed=True, message="executed")


def test_unknown_requirement_cannot_inherit_json_success():
    failures = missing_required_checks(["must_run_custom_check"], [check("json_valid_payload")])
    assert len(failures) == 1 and not failures[0].passed


def test_alias_requires_the_actual_equivalent_check():
    assert not missing_required_checks(["finite_expression_check"], [check("expression_finite_values")])
    assert missing_required_checks(["finite_expression_check"], [check("table_finite_values")])


def test_aggregate_auditor_rejects_unexecuted_required_check(tmp_path):
    from eacbp.auditor import ScientificAuditor
    from eacbp.artifact.registry import ArtifactRegistry
    from eacbp.schemas.artifact import ArtifactType
    from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus
    registry = ArtifactRegistry(str(tmp_path))
    uri = "json://s/result/v1"
    registry.register(uri, {"valid": True}, ArtifactType.JSON, "s", "custom", "custom")
    task = TaskContract(task_id="custom", capability="custom", validation_requirements=["must_run_custom_check"])
    result = TaskResult(task_id="custom", capability="custom", method_used="custom", status=TaskStatus.SUCCESS, output_artifacts=[uri])
    report = ScientificAuditor().audit_task(task, result, registry)
    assert not report.overall_passed
    assert any(c.check_name == "json_valid_payload" and c.passed for c in report.checks)
    assert any(c.check_name == "required_check_missing:must_run_custom_check" and not c.passed for c in report.checks)
