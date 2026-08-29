"""
Side Effect Validator ensuring that capability execution strictly obeys the TaskContract.
Prevents agents or scripts from silently re-clustering, re-filtering, or mutating parent datasets.
"""

from typing import Tuple, List, Dict, Any, Optional
from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus, ExecutionFailureType
from eacbp.artifact.registry import ArtifactRegistry


class SideEffectValidator:
    """Validates that a TaskResult strictly adheres to allowed/forbidden operations and bounds."""

    @staticmethod
    def validate(
        contract: TaskContract,
        result: TaskResult,
        input_payloads: Optional[Dict[str, Any]] = None,
        output_payloads: Optional[Dict[str, Any]] = None,
    ) -> Tuple[bool, Optional[str]]:
        # 1. Check operation list violation
        for op in result.executed_operations:
            if contract.forbidden_operations and op in contract.forbidden_operations:
                return False, f"POLICY VIOLATION: Operation '{op}' is explicitly FORBIDDEN in TaskContract {contract.task_id}."

        if contract.allowed_operations:
            unallowed = [op for op in result.executed_operations if op not in contract.allowed_operations]
            if unallowed:
                return False, f"POLICY VIOLATION: Operations {unallowed} were executed but not in allowed_operations list."

        # 2. Check data invariant side effects
        if input_payloads and output_payloads and contract.input_artifacts:
            in_key = contract.input_artifacts[0]
            in_data = input_payloads.get(in_key)
            if in_data is not None and result.output_artifacts:
                out_key = result.output_artifacts[0]
                out_data = output_payloads.get(out_key)
                
                # Check cell count mutation if filter_cells is forbidden
                if "filter_cells" in contract.forbidden_operations:
                    in_n_cells = getattr(in_data, "n_obs", None)
                    if in_n_cells is None and isinstance(in_data, dict) and ("obs" in in_data or "X" in in_data):
                        in_n_cells = len(in_data.get("obs", [])) if "obs" in in_data else getattr(in_data.get("X"), "shape", [0])[0]

                    out_n_cells = getattr(out_data, "n_obs", None)
                    if out_n_cells is None and isinstance(out_data, dict) and ("obs" in out_data or "X" in out_data):
                        out_n_cells = len(out_data.get("obs", [])) if "obs" in out_data else getattr(out_data.get("X"), "shape", [0])[0]

                    if in_n_cells is not None and in_n_cells > 0 and out_n_cells is not None and in_n_cells != out_n_cells:
                        return False, f"POLICY VIOLATION: Cell count changed from {in_n_cells} to {out_n_cells} despite 'filter_cells' being forbidden."

                # Check cluster label tampering if recluster is forbidden
                if "recluster" in contract.forbidden_operations and hasattr(in_data, "obs") and hasattr(out_data, "obs"):
                    if "leiden" in in_data.obs.columns and "leiden" in out_data.obs.columns:
                        if not in_data.obs["leiden"].equals(out_data.obs["leiden"]):
                            return False, f"POLICY VIOLATION: Cluster assignments 'leiden' were modified despite 'recluster' being forbidden."

        return True, None
