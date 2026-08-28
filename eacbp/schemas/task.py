"""
Task Contract and Result schemas for controlling capability execution and preventing unauthorized upstream modifications.
"""

from typing import List, Dict, Any, Optional
from enum import Enum
from pydantic import BaseModel, Field


class TaskStatus(str, Enum):
    PENDING = "pending"
    RUNNING = "running"
    SUCCESS = "success"
    EXECUTION_FAILURE = "execution_failure"
    METHOD_FAILURE = "method_failure"
    SCIENTIFIC_FAILURE = "scientific_failure"
    POLICY_VIOLATION = "policy_violation"


class ExecutionFailureType(str, Enum):
    CODE_ERROR = "code_error"
    MEMORY_ERROR = "memory_error"
    DEPENDENCY_ERROR = "dependency_error"
    CONVERGENCE_ERROR = "convergence_error"
    INSTABILITY_ERROR = "instability_error"
    INSUFFICIENT_EVIDENCE = "insufficient_evidence"
    UNAUTHORIZED_SIDE_EFFECT = "unauthorized_side_effect"


class RetryPolicy(BaseModel):
    max_execution_retry: int = Field(2, description="Max retries for execution / runtime failures")
    max_method_retry: int = Field(2, description="Max retries for method-level algorithm fallback")
    fallback_methods: List[str] = Field(default_factory=list, description="Ordered list of alternative fallback methods")
    require_human_after: int = Field(4, description="Escalate to human review after N cumulative attempts")


class TaskContract(BaseModel):
    task_id: str = Field(..., description="Unique task identifier, e.g., task_018")
    capability: str = Field(..., description="Target capability name, e.g., trajectory_inference")
    method: Optional[str] = Field(None, description="Requested method implementation, e.g., cellrank, paga")
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
    execution_time_sec: float = Field(0.0)
    logs: str = Field("", description="Captured stdout / execution log stream")
    error_type: Optional[ExecutionFailureType] = Field(None)
    error_message: Optional[str] = Field(None)
