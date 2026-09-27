"""
Task Contract and Result schemas for controlling capability execution and preventing unauthorized upstream modifications.
"""

from typing import List, Dict, Any, Optional
from enum import Enum
import hashlib
import json

from pydantic import BaseModel, Field, model_validator


class TaskStatus(str, Enum):
    PENDING = "pending"
    RUNNING = "running"
    SUCCESS = "success"
    EXECUTION_FAILURE = "execution_failure"
    METHOD_FAILURE = "method_failure"
    SCIENTIFIC_FAILURE = "scientific_failure"
    POLICY_VIOLATION = "policy_violation"
    BLOCKED = "blocked"


class ExecutionFailureType(str, Enum):
    TIMEOUT = "timeout"
    CANCELLED = "cancelled"
    RESOURCE_POLICY = "resource_policy"
    CODE_ERROR = "code_error"
    MEMORY_ERROR = "memory_error"
    DEPENDENCY_ERROR = "dependency_error"
    CONVERGENCE_ERROR = "convergence_error"
    INSTABILITY_ERROR = "instability_error"
    INSUFFICIENT_EVIDENCE = "insufficient_evidence"
    UNAUTHORIZED_SIDE_EFFECT = "unauthorized_side_effect"


class ContrastSpec(BaseModel):
    """Canonical resolved contrast shared by statistical contracts and outputs."""

    kind: str = Field(..., description="categorical or numeric contrast")
    vector: tuple[Any, ...] = Field(..., description="Contrast passed to the statistical method")
    factor: Optional[str] = None
    tested_level: Optional[str] = None
    reference_level: Optional[str] = None
    label: str
    effect_definition: str

    def to_pydeseq2(self):
        if self.kind == "numeric":
            import numpy as np

            return np.asarray(self.vector, dtype=float)
        return [str(value) for value in self.vector]

    def as_dict(self) -> Dict[str, Any]:
        return {
            "kind": self.kind,
            "vector": list(self.vector),
            "factor": self.factor,
            "tested_level": self.tested_level,
            "reference_level": self.reference_level,
            "label": self.label,
            "effect_definition": self.effect_definition,
        }


class AssumptionStatus(str, Enum):
    PASSED = "passed"
    FAILED = "failed"
    UNKNOWN = "unknown"
    NOT_ASSESSABLE = "not_assessable"


class AssumptionAssessment(BaseModel):
    name: str
    status: AssumptionStatus
    scope: str = "analysis"
    evidence: Dict[str, Any] = Field(default_factory=dict)
    reason: str = ""


class ScientificResultStatus(str, Enum):
    ESTIMATED_SUPPORTED = "estimated_supported"
    ESTIMATED_INCONCLUSIVE = "estimated_inconclusive"
    NOT_ESTIMABLE = "not_estimable"
    ASSUMPTIONS_FAILED = "assumptions_failed"


class ScientificResult(BaseModel):
    """Neutral account of whether a requested statistical result was estimated."""

    status: ScientificResultStatus
    summary: str
    n_features_tested: int = Field(0, ge=0)
    n_features_supported: int = Field(0, ge=0)
    assumptions: List[AssumptionAssessment] = Field(default_factory=list)

    @model_validator(mode="after")
    def validate_counts(self):
        if self.n_features_supported > self.n_features_tested:
            raise ValueError("supported features cannot exceed tested features")
        return self


class InferenceContract(BaseModel):
    """Versioned, fingerprinted description of the requested estimand and design."""

    contract_version: int = Field(1, ge=1)
    contract_id: str = ""
    scientific_question: str
    estimand: str
    target_population: str
    target_cell_type: Optional[str] = None
    feature_scope: str
    input_artifact_uris: List[str] = Field(default_factory=list)
    independent_unit: str
    observation_unit: str
    paired: bool
    design_formula: str
    counts_source: str
    contrast_spec: ContrastSpec
    alpha: float = Field(0.05, gt=0, lt=1)
    confidence_level: float = Field(0.95, gt=0, lt=1)
    fdr_method: str = "benjamini_hochberg"
    fdr_family: str
    sensitivity_plan: List[str] = Field(default_factory=list)

    @model_validator(mode="after")
    def validate_identity_and_interval(self):
        if abs(self.confidence_level - (1.0 - self.alpha)) > 1e-12:
            raise ValueError("confidence_level must equal 1 - alpha")
        payload = self.model_dump(mode="json", exclude={"contract_id"})
        canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False)
        expected = f"inference-v{self.contract_version}:{hashlib.sha256(canonical.encode('utf-8')).hexdigest()}"
        if self.contract_id and self.contract_id != expected:
            raise ValueError("inference contract_id does not match its versioned content")
        self.contract_id = expected
        return self


class RetryPolicy(BaseModel):
    max_execution_retry: int = Field(2, ge=0, description="Max retries for execution / runtime failures")
    max_method_retry: int = Field(2, ge=0, description="Max retries for method-level algorithm fallback")
    fallback_methods: List[str] = Field(default_factory=list, description="Ordered list of alternative fallback methods")
    require_human_after: int = Field(4, ge=1, description="Maximum attempts before recording a failure requiring review")


class TaskContract(BaseModel):
    task_id: str = Field(..., description="Unique task identifier, e.g., task_018")
    capability: str = Field(..., description="Target capability name, e.g., trajectory_inference")
    method: Optional[str] = Field(None, description="Requested method implementation, e.g., cellrank, paga")
    depends_on: List[str] = Field(default_factory=list)
    input_artifacts: List[str] = Field(default_factory=list, description="Input artifact URIs, e.g., ['adata://AD/microglia/v4']")
    
    # Contract bounds preventing agent rogue upstream alterations
    allowed_operations: List[str] = Field(
        default_factory=list,
        description="Explicitly allowed operations, e.g., ['build_neighbor_graph', 'infer_trajectory']"
    )
    forbidden_operations: List[str] = Field(
        default_factory=list,
        description="Explicitly forbidden operations, e.g., ['filter_cells', 'normalize', 'batch_correct', 'recluster']"
    )
    
    parameters: Dict[str, Any] = Field(default_factory=dict, description="Execution parameters and hyperparams")
    expected_outputs: List[str] = Field(default_factory=list, description="Expected output artifact categories")
    validation_requirements: List[str] = Field(
        default_factory=list,
        description="Validation checks required, e.g. ['topology_stability', 'root_sensitivity', 'marker_consistency']"
    )
    inference_contract: Optional[InferenceContract] = Field(
        None, description="Resolved scientific estimand/design contract, included in resume signatures"
    )
    retry_policy: RetryPolicy = Field(default_factory=RetryPolicy)


class TaskResult(BaseModel):
    task_id: str = Field(...)
    status: TaskStatus = Field(...)
    capability: str = Field(...)
    method_used: str = Field(...)
    input_artifacts: List[str] = Field(default_factory=list)
    output_artifacts: List[str] = Field(default_factory=list, description="Produced artifact URIs")
    executed_operations: List[str] = Field(default_factory=list, description="List of operations actually performed")
    metrics: Dict[str, Any] = Field(default_factory=dict, description="Summary quantitative metrics from run")
    inference_contract: Optional[InferenceContract] = None
    inference_contract_id: Optional[str] = None
    scientific_result: Optional[ScientificResult] = None
    execution_time_sec: float = Field(0.0)
    logs: str = Field("", description="Captured stdout / execution log stream")
    error_type: Optional[ExecutionFailureType] = Field(None)
    error_message: Optional[str] = Field(None)

    @model_validator(mode="after")
    def validate_inference_contract_identity(self):
        if self.inference_contract is not None:
            contract_id = self.inference_contract.contract_id
            if self.inference_contract_id is not None and self.inference_contract_id != contract_id:
                raise ValueError("TaskResult inference_contract_id does not match inference_contract")
            self.inference_contract_id = contract_id
        return self
