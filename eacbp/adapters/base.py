"""
Base Agent Adapter interface and contract guardrails integration.
Enforces TaskContract bounds, sandboxing, pre-execution snapshots, and versioned artifact wrapping.
"""

from abc import ABC, abstractmethod
from typing import Dict, List, Any, Optional, Tuple
import time
import copy
import hashlib
import json
import numpy as np
import pandas as pd

from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus, ExecutionFailureType
from eacbp.schemas.artifact import ArtifactType, ArtifactMetadata
from eacbp.capabilities.base import BaseCapability, ImplementationType
from eacbp.capabilities.side_effect import SideEffectValidator
from eacbp.capabilities.sc_data import SCData
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.artifact.uri import ArtifactURI


class BaseAgentAdapter(BaseCapability, ABC):
    """
    Abstract base adapter for third-party biological agents.
    Integrates seamlessly with TaskContract and SideEffectValidator.
    
    Provides:
    1. Pre-execution input validation and immutable state snapshotting.
    2. Sandboxed execution dispatching to concrete _execute_agent logic.
    3. Strict post-execution contract guardrail validation (intercepting rogue side effects).
    4. Versioned, SHA-256 hashed artifact wrapping with complete provenance lineage.
    """

    def __init__(
        self,
        capability_name: str,
        implementation_id: str,
        accepts_modalities: Optional[List[str]] = None,
        accepts_types: Optional[List[ArtifactType]] = None,
        requires_keys: Optional[List[str]] = None,
        suitable_for: Optional[List[str]] = None,
        output_types: Optional[List[ArtifactType]] = None,
        agent_config: Optional[Dict[str, Any]] = None,
    ):
        super().__init__(
            capability_name=capability_name,
            implementation_id=implementation_id,
            implementation_type=ImplementationType.AGENT_ADAPTER,
            accepts_modalities=accepts_modalities or ["scRNA"],
            accepts_types=accepts_types or [ArtifactType.ANNDATA],
            requires_keys=requires_keys or [],
            suitable_for=suitable_for or [],
            output_types=output_types or [ArtifactType.TABLE, ArtifactType.JSON],
        )
        self.agent_config = agent_config or {}

    def execute(self, contract: TaskContract, registry: ArtifactRegistry) -> TaskResult:
        """
        Executes the agent adapter under strict TaskContract bounds.
        Performs pre-execution snapshotting, agent invocation, and post-execution guardrail verification.
        """
        start_time = time.time()

        # 1. Verify and load input artifacts
        input_payloads: Dict[str, Any] = {}
        for in_uri in contract.input_artifacts:
            if not registry.exists(in_uri):
                return TaskResult(
                    task_id=contract.task_id,
                    status=TaskStatus.EXECUTION_FAILURE,
                    capability=self.capability_name,
                    method_used=self.implementation_id,
                    input_artifacts=contract.input_artifacts,
                    output_artifacts=[],
                    executed_operations=[],
                    execution_time_sec=time.time() - start_time,
                    error_type=ExecutionFailureType.DEPENDENCY_ERROR,
                    error_message=f"Input artifact '{in_uri}' not found in ArtifactRegistry.",
                )
            input_payloads[in_uri] = registry.load_payload(in_uri)

        # 2. Capture pre-execution fingerprint of input payloads to detect in-place mutations
        pre_execution_fingerprints = {
            uri: self._compute_payload_fingerprint(payload)
            for uri, payload in input_payloads.items()
        }

        # 3. Execute concrete agent implementation
        try:
            task_result = self._execute_agent(contract, registry, input_payloads)
        except Exception as e:
            return TaskResult(
                task_id=contract.task_id,
                status=TaskStatus.EXECUTION_FAILURE,
                capability=self.capability_name,
                method_used=self.implementation_id,
                input_artifacts=contract.input_artifacts,
                output_artifacts=[],
                executed_operations=[],
                execution_time_sec=time.time() - start_time,
                error_type=ExecutionFailureType.CODE_ERROR,
                error_message=f"Adapter execution crashed with exception: {str(e)}",
            )

        task_result.execution_time_sec = time.time() - start_time

        # 4. Invariant & Contract Guardrail Validations
        violation_reason = self._check_guardrail_violations(
            contract=contract,
            result=task_result,
            input_payloads=input_payloads,
            pre_fingerprints=pre_execution_fingerprints,
            registry=registry,
        )

        if violation_reason is not None:
            task_result.status = TaskStatus.POLICY_VIOLATION
            task_result.error_type = ExecutionFailureType.UNAUTHORIZED_SIDE_EFFECT
            task_result.error_message = violation_reason

        return task_result

    @abstractmethod
    def _execute_agent(
        self,
        contract: TaskContract,
        registry: ArtifactRegistry,
        input_payloads: Dict[str, Any],
    ) -> TaskResult:
        """
        Concrete adapter logic. Must be implemented by specialized agent adapters.
        """
        pass

    @staticmethod
    def _hash_ndarray(arr: np.ndarray, hasher: Any) -> None:
        """Deterministically hashes a numpy ndarray."""
        hasher.update(b"__NDARRAY__")
        hasher.update(str(arr.shape).encode("utf-8"))
        hasher.update(str(arr.dtype).encode("utf-8"))
        if arr.size > 0:
            if not (arr.flags.c_contiguous or arr.flags.f_contiguous):
                arr = np.ascontiguousarray(arr)
            hasher.update(arr.tobytes())

    @staticmethod
    def _hash_dataframe(df: pd.DataFrame, hasher: Any) -> None:
        """Deterministically hashes a pandas DataFrame including columns, index, dtypes, and values."""
        hasher.update(b"__DATAFRAME__")
        hasher.update(str(df.shape).encode("utf-8"))
        hasher.update(str(list(df.columns)).encode("utf-8"))
        hasher.update(str(list(df.index)).encode("utf-8"))
        hasher.update(str(df.dtypes.to_dict()).encode("utf-8"))
        if not df.empty:
            try:
                hasher.update(pd.util.hash_pandas_object(df, index=True).values.tobytes())
            except Exception:
                try:
                    hasher.update(df.to_json(orient="split", default_handler=str).encode("utf-8"))
                except Exception:
                    hasher.update(df.to_csv().encode("utf-8"))

    @staticmethod
    def _hash_series(s: pd.Series, hasher: Any) -> None:
        """Deterministically hashes a pandas Series."""
        hasher.update(b"__SERIES__")
        hasher.update(str(s.name).encode("utf-8"))
        hasher.update(str(s.shape).encode("utf-8"))
        hasher.update(str(s.dtype).encode("utf-8"))
        if not s.empty:
            try:
                hasher.update(pd.util.hash_pandas_object(s, index=True).values.tobytes())
            except Exception:
                hasher.update(str(s.to_dict()).encode("utf-8"))

    def _hash_recursive(self, obj: Any, hasher: Any, depth: int = 0) -> None:
        """Recursively hashes an arbitrary payload object."""
        if depth > 50:
            hasher.update(str(obj).encode("utf-8"))
            return

        if obj is None:
            hasher.update(b"__NONE__")
        elif isinstance(obj, SCData):
            hasher.update(b"__SCDATA__")
            self._hash_ndarray(obj.X, hasher)
            self._hash_dataframe(obj.obs, hasher)
            self._hash_dataframe(obj.var, hasher)
            if hasattr(obj, "obsm") and obj.obsm:
                hasher.update(b"__OBSM__")
                for k in sorted(obj.obsm.keys(), key=lambda x: str(x)):
                    hasher.update(str(k).encode("utf-8"))
                    self._hash_recursive(obj.obsm[k], hasher, depth + 1)
            if hasattr(obj, "obsp") and obj.obsp:
                hasher.update(b"__OBSP__")
                for k in sorted(obj.obsp.keys(), key=lambda x: str(x)):
                    hasher.update(str(k).encode("utf-8"))
                    self._hash_recursive(obj.obsp[k], hasher, depth + 1)
            if hasattr(obj, "uns") and obj.uns:
                hasher.update(b"__UNS__")
                for k in sorted(obj.uns.keys(), key=lambda x: str(x)):
                    hasher.update(str(k).encode("utf-8"))
                    self._hash_recursive(obj.uns[k], hasher, depth + 1)
        elif isinstance(obj, pd.DataFrame):
            self._hash_dataframe(obj, hasher)
        elif isinstance(obj, pd.Series):
            self._hash_series(obj, hasher)
        elif isinstance(obj, np.ndarray):
            self._hash_ndarray(obj, hasher)
        elif isinstance(obj, dict):
            hasher.update(b"__DICT__")
            for k in sorted(obj.keys(), key=lambda x: str(x)):
                hasher.update(str(k).encode("utf-8"))
                self._hash_recursive(obj[k], hasher, depth + 1)
        elif isinstance(obj, (list, tuple)):
            hasher.update(f"__{type(obj).__name__.upper()}__".encode("utf-8"))
            hasher.update(str(len(obj)).encode("utf-8"))
            for item in obj:
                self._hash_recursive(item, hasher, depth + 1)
        elif isinstance(obj, (set, frozenset)):
            hasher.update(b"__SET__")
            hasher.update(str(len(obj)).encode("utf-8"))
            for item in sorted(list(obj), key=lambda x: str(x)):
                self._hash_recursive(item, hasher, depth + 1)
        elif isinstance(obj, (bytes, bytearray)):
            hasher.update(b"__BYTES__")
            hasher.update(bytes(obj))
        elif isinstance(obj, (int, float, str, bool)):
            hasher.update(f"__{type(obj).__name__.upper()}__{obj}".encode("utf-8"))
        else:
            hasher.update(str(obj).encode("utf-8"))

    def _compute_payload_fingerprint(self, payload: Any) -> str:
        """Computes a comprehensive SHA-256 fingerprint of payload state to detect in-place mutation."""
        hasher = hashlib.sha256()
        self._hash_recursive(payload, hasher)
        return hasher.hexdigest()

    def _check_guardrail_violations(
        self,
        contract: TaskContract,
        result: TaskResult,
        input_payloads: Dict[str, Any],
        pre_fingerprints: Dict[str, str],
        registry: ArtifactRegistry,
    ) -> Optional[str]:
        """
        Comprehensive contract guardrail checks.
        Returns error message if violation is detected, None if valid.
        """
        # A. Check executed operations vs forbidden_operations
        for op in result.executed_operations:
            if contract.forbidden_operations and op in contract.forbidden_operations:
                return (
                    f"POLICY VIOLATION: Operation '{op}' is explicitly FORBIDDEN in "
                    f"TaskContract {contract.task_id}."
                )

        # B. Check executed operations vs allowed_operations
        if contract.allowed_operations:
            unallowed = [op for op in result.executed_operations if op not in contract.allowed_operations]
            if unallowed:
                return (
                    f"POLICY VIOLATION: Operations {unallowed} were executed but not in "
                    f"allowed_operations list for TaskContract {contract.task_id}."
                )

        # C. Check in-place mutation of input payloads
        if "in_place_mutation" in contract.forbidden_operations or "mutate_raw" in contract.forbidden_operations:
            for uri, in_data in input_payloads.items():
                post_fp = self._compute_payload_fingerprint(in_data)
                pre_fp = pre_fingerprints.get(uri)
                if pre_fp and post_fp != pre_fp:
                    return (
                        f"POLICY VIOLATION: Input artifact '{uri}' was mutated in-place "
                        f"despite in-place mutation being forbidden in TaskContract {contract.task_id}."
                    )

        # D. Check output payloads for silent cell filtering and cluster tampering
        output_payloads: Dict[str, Any] = {}
        for out_uri in result.output_artifacts:
            if registry.exists(out_uri):
                output_payloads[out_uri] = registry.load_payload(out_uri)

        # Check data invariants (handling both SCData and dict formats)
        if input_payloads and output_payloads and contract.input_artifacts:
            in_key = contract.input_artifacts[0]
            in_data = input_payloads.get(in_key)
            if in_data is not None and result.output_artifacts:
                out_key = result.output_artifacts[0]
                out_data = output_payloads.get(out_key)

                in_obs = in_data.obs if hasattr(in_data, "obs") else (in_data.get("obs") if isinstance(in_data, dict) else None)
                out_obs = out_data.obs if hasattr(out_data, "obs") else (out_data.get("obs") if isinstance(out_data, dict) else None)

                # Check cell count mutation if filter_cells is forbidden
                if "filter_cells" in contract.forbidden_operations:
                    in_n_cells = getattr(in_data, "n_obs", None) or (len(in_obs) if in_obs is not None else None)
                    out_n_cells = getattr(out_data, "n_obs", None) or (len(out_obs) if out_obs is not None else None)
                    if in_n_cells is not None and out_n_cells is not None and in_n_cells != out_n_cells:
                        return f"POLICY VIOLATION: Cell count changed from {in_n_cells} to {out_n_cells} despite 'filter_cells' being forbidden."

                # Check cluster label tampering if recluster is forbidden
                if "recluster" in contract.forbidden_operations or "mutate_clusters" in contract.forbidden_operations:
                    if in_obs is not None and out_obs is not None:
                        for cluster_col in ["leiden", "louvain", "cluster", "cell_type", "cell_type_annotation"]:
                            if cluster_col in in_obs.columns and cluster_col in out_obs.columns:
                                if not in_obs[cluster_col].equals(out_obs[cluster_col]):
                                    return f"POLICY VIOLATION: Cluster assignments '{cluster_col}' were modified despite 'recluster' being forbidden."

        # Run standard SideEffectValidator check
        valid, violation_msg = SideEffectValidator.validate(
            contract=contract,
            result=result,
            input_payloads=input_payloads,
            output_payloads=output_payloads,
        )
        if not valid:
            return violation_msg

        return None

    def _to_sc_data(self, payload: Any) -> SCData:
        """Utility to convert dict / SCData payload into SCData instance."""
        if isinstance(payload, SCData):
            return payload
        if isinstance(payload, dict):
            return SCData.from_dict(payload)
        raise TypeError(f"Cannot convert payload of type {type(payload)} to SCData.")

    def _generate_output_uri(
        self,
        study_id: str,
        stage: str,
        scheme: str = "adata",
        version: str = "v1",
    ) -> str:
        """Generates a canonical Artifact URI string."""
        return f"{scheme}://{study_id}/{stage}/{version}"

    def _register_versioned_artifact(
        self,
        registry: ArtifactRegistry,
        uri_str: str,
        payload: Any,
        artifact_type: ArtifactType,
        study_id: str,
        task_id: str,
        operation: str,
        parent_uris: Optional[List[str]] = None,
        parameters: Optional[Dict[str, Any]] = None,
        summary_metrics: Optional[Dict[str, Any]] = None,
    ) -> ArtifactMetadata:
        """Registers an artifact with full lineage tracking and SHA-256 content addressing."""
        return registry.register(
            uri_str=uri_str,
            payload=payload,
            artifact_type=artifact_type,
            study_id=study_id,
            created_by_task=task_id,
            operation=operation,
            parent_uris=parent_uris or [],
            parameters=parameters or {},
            software_versions={
                "adapter": self.implementation_id,
                "capability": self.capability_name,
                "eacbp": "0.1.0",
            },
            summary_metrics=summary_metrics or {},
        )
