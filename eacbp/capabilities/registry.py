"""
Global Capability Registry for discovering and resolving computational implementations.
"""

from typing import Dict, List, Optional, Any
from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus, ExecutionFailureType
from eacbp.capabilities.base import BaseCapability, ImplementationType
from eacbp.capabilities.side_effect import SideEffectValidator
from eacbp.artifact.registry import ArtifactRegistry


class CapabilityRegistry:
    """Registry maintaining all available computational capabilities across engines and adapters."""

    def __init__(self):
        # Key: (capability_name, implementation_id) -> BaseCapability
        self._capabilities: Dict[str, Dict[str, BaseCapability]] = {}

    def register(self, capability: BaseCapability):
        cap_name = capability.capability_name
        impl_id = capability.implementation_id
        if cap_name not in self._capabilities:
            self._capabilities[cap_name] = {}
        self._capabilities[cap_name][impl_id] = capability

    def get(self, capability_name: str, method: Optional[str] = None) -> BaseCapability:
        if capability_name not in self._capabilities:
            raise KeyError(f"No implementations registered for capability: '{capability_name}'")
        
        impls = self._capabilities[capability_name]
        if method:
            if method not in impls:
                raise KeyError(f"Method '{method}' not found for capability '{capability_name}'. Available: {list(impls.keys())}")
            return impls[method]
        
        # Default to the first registered implementation
        return next(iter(impls.values()))

    def list_capabilities(self) -> Dict[str, List[Dict[str, Any]]]:
        result = {}
        for cap_name, impls in self._capabilities.items():
            result[cap_name] = [
                {
                    "implementation_id": cap.implementation_id,
                    "type": cap.implementation_type.value,
                    "suitable_for": cap.suitable_for,
                    "accepts_modalities": cap.accepts_modalities,
                }
                for cap in impls.values()
            ]
        return result

    def execute_contract(self, contract: TaskContract, registry: ArtifactRegistry) -> TaskResult:
        """Executes a TaskContract through the registered capability with side-effect validation."""
        capability = self.get(contract.capability, contract.method)

        # Load input payloads for side-effect pre-check
        input_payloads = {}
        for in_uri in contract.input_artifacts:
            if registry.exists(in_uri):
                input_payloads[in_uri] = registry.load_payload(in_uri)

        # Run capability execution
        task_result = capability.execute(contract, registry)

        # If execution succeeded, run side-effect validation against contract
        if task_result.status == TaskStatus.SUCCESS:
            output_payloads = {}
            for out_uri in task_result.output_artifacts:
                if registry.exists(out_uri):
                    output_payloads[out_uri] = registry.load_payload(out_uri)

            valid, violation_msg = SideEffectValidator.validate(
                contract=contract,
                result=task_result,
                input_payloads=input_payloads,
                output_payloads=output_payloads,
            )

            if not valid:
                task_result.status = TaskStatus.POLICY_VIOLATION
                task_result.error_type = ExecutionFailureType.UNAUTHORIZED_SIDE_EFFECT
                task_result.error_message = violation_msg

        return task_result
