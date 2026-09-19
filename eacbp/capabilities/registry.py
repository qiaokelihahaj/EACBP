"""
Global Capability Registry for discovering and resolving computational implementations.
"""

from typing import Dict, List, Optional, Any, Iterable
from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus, ExecutionFailureType
from eacbp.capabilities.base import BaseCapability, CapabilityDescriptor, ImplementationType
from eacbp.capabilities.side_effect import SideEffectValidator
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.schemas.artifact import ArtifactType
from eacbp.artifact.uri import ArtifactURI


class CapabilityRegistry:
    """Registry maintaining all available computational capabilities across engines and adapters."""

    def __init__(self):
        # Key: (capability_name, implementation_id) -> BaseCapability
        self._capabilities: Dict[str, Dict[str, BaseCapability]] = {}

    def register(self, capability: BaseCapability, *, overwrite: bool = False):
        """Register one capability implementation.

        Capability IDs are part of the execution contract.  Accidentally
        replacing an implementation under the same ``(capability, method)``
        key can therefore make a caller run a different algorithm than the
        one it requested.  Registration is deliberately strict by default;
        callers that really intend to replace an implementation must say so
        with ``overwrite=True``.

        The method returns the registered capability so factories and plugin
        loaders can compose registration without reaching into the private
        mapping.
        """
        if not isinstance(capability, BaseCapability):
            raise TypeError("CapabilityRegistry.register expects a BaseCapability instance")
        descriptor = capability.get_descriptor()
        if descriptor.capability_name != capability.capability_name:
            raise ValueError(
                "Capability descriptor name does not match implementation: "
                f"{descriptor.capability_name!r} != {capability.capability_name!r}"
            )
        if descriptor.method != capability.implementation_id:
            raise ValueError(
                "Capability descriptor method does not match implementation: "
                f"{descriptor.method!r} != {capability.implementation_id!r}"
            )
        cap_name = capability.capability_name
        impl_id = capability.implementation_id
        if cap_name not in self._capabilities:
            self._capabilities[cap_name] = {}
        if impl_id in self._capabilities[cap_name] and not overwrite:
            raise ValueError(
                f"Capability implementation already registered: "
                f"'{cap_name}/{impl_id}'. Pass overwrite=True to replace it."
            )
        self._capabilities[cap_name][impl_id] = capability
        return capability

    def has(self, capability_name: str, method: Optional[str] = None) -> bool:
        """Return whether a capability (or a specific implementation) exists.

        Alias resolution is intentionally delegated to :meth:`get`; this
        helper is useful for assembly code that only needs to test presence.
        """
        try:
            self.get(capability_name, method)
        except KeyError:
            return False
        return True

    def get(self, capability_name: str, method: Optional[str] = None) -> BaseCapability:
        if capability_name not in self._capabilities:
            raise KeyError(f"No implementations registered for capability: '{capability_name}'")
        
        impls = self._capabilities[capability_name]
        if method:
            if method not in impls:
                for capability in impls.values():
                    descriptor = capability.get_descriptor()
                    aliases = set(getattr(capability, "legacy_aliases", {})) | set(descriptor.method_aliases)
                    if method in aliases:
                        return capability
                raise KeyError(f"Method '{method}' not found for capability '{capability_name}'. Available: {list(impls.keys())}")
            return impls[method]
        
        # Default to the first registered implementation
        return next(iter(impls.values()))

    def resolve(self, capability_name: str, method: Optional[str] = None) -> BaseCapability:
        """Resolve a capability implementation, including declared aliases."""

        return self.get(capability_name, method)

    def describe(self, capability_name: str, method: Optional[str] = None) -> CapabilityDescriptor:
        """Return the executable descriptor for a registered implementation."""

        return self.resolve(capability_name, method).get_descriptor()

    @staticmethod
    def _artifact_type_from_uri(uri: str) -> Optional[ArtifactType]:
        """Map canonical artifact URI schemes to declared artifact types."""

        try:
            scheme = ArtifactURI.parse(uri).scheme
        except (TypeError, ValueError):
            return None
        schemes = {
            "adata": ArtifactType.ANNDATA,
            "anndata": ArtifactType.ANNDATA,
            "spatial": ArtifactType.SPATIAL_DATA,
            "table": ArtifactType.TABLE,
            "json": ArtifactType.JSON,
            "fastq": ArtifactType.FASTQ,
            "fig": ArtifactType.FIGURE,
            "figure": ArtifactType.FIGURE,
            "graph": ArtifactType.GRAPH,
            "report": ArtifactType.REPORT,
            "genes": ArtifactType.GENE_LIST,
            "gene_list": ArtifactType.GENE_LIST,
        }
        return schemes.get(scheme)

    def validate_contract(
        self,
        contract: TaskContract,
        artifact_registry: Optional[ArtifactRegistry] = None,
        *,
        add_required_audits: bool = True,
    ) -> TaskContract:
        """Validate and normalise a task contract against its descriptor.

        ``input_types`` and ``output_types`` are alternatives: an input is
        valid when its metadata type is one of the declared types.  Output
        URI schemes are checked statically; output payload metadata is checked
        by the executor after the capability has run.  The method mutates and
        returns the supplied Pydantic contract for convenient planner use.
        """

        if not isinstance(contract, TaskContract):
            raise TypeError("validate_contract expects a TaskContract")
        capability = self.resolve(contract.capability, contract.method)
        descriptor = capability.get_descriptor()
        contract.parameters = descriptor.validate_parameters(contract.parameters)
        if contract.method is None:
            contract.method = capability.implementation_id
        if add_required_audits:
            contract.validation_requirements = list(dict.fromkeys(
                [*contract.validation_requirements, *descriptor.required_audit_ids]
            ))

        if descriptor.validate_types and descriptor.output_types:
            invalid_outputs = []
            for uri in contract.expected_outputs:
                artifact_type = self._artifact_type_from_uri(uri)
                if artifact_type is not None and artifact_type not in descriptor.output_types:
                    invalid_outputs.append((uri, artifact_type.value))
            if invalid_outputs:
                allowed = [value.value for value in descriptor.output_types]
                raise ValueError(
                    f"Capability '{descriptor.capability_name}/{descriptor.method}' declares output types "
                    f"{allowed}, received incompatible outputs: {invalid_outputs}"
                )

        if descriptor.validate_types and artifact_registry is not None and descriptor.input_types:
            missing = [uri for uri in contract.input_artifacts if not artifact_registry.exists(uri)]
            if missing:
                raise KeyError(f"Capability input artifacts are not registered: {missing}")
            incompatible = []
            for uri in contract.input_artifacts:
                metadata = artifact_registry.get_metadata(uri)
                if metadata.type not in descriptor.input_types:
                    incompatible.append((uri, metadata.type.value))
            if incompatible:
                allowed = [value.value for value in descriptor.input_types]
                raise ValueError(
                    f"Capability '{descriptor.capability_name}/{descriptor.method}' accepts input types "
                    f"{allowed}, received incompatible inputs: {incompatible}"
                )
        return contract

    # ``prepare_contract`` reads naturally at execution boundaries and keeps
    # plugin code independent of the registry's historical method name.
    prepare_contract = validate_contract

    def audit_validators(self):
        from eacbp.auditor.descriptor import DescriptorValidator
        return [DescriptorValidator(cap.get_descriptor())
                for methods in self._capabilities.values() for cap in methods.values()
                if cap.get_descriptor().validator is not None]

    def plan_extensions(self, manifest, state, tasks):
        """Build user-requested extension tasks from registered declarations."""
        from eacbp.capabilities.advanced_descriptors import ADVANCED_EXTENSION_DESCRIPTORS
        built = []
        for name, settings in state.get("analysis_extensions", {}).items():
            if name in ADVANCED_EXTENSION_DESCRIPTORS or settings is False:
                continue
            if settings is not True and not isinstance(settings, dict):
                raise ValueError(f"Extension {name} must be true, false, or a parameter object")
            descriptor = self.describe(name, state.get("method_overrides", {}).get(name))
            if descriptor.plan_factory is None:
                raise ValueError(f"Extension {name!r} has no declared plan factory")
            values = dict(settings) if isinstance(settings, dict) else {}
            values.update(state.get("capability_parameters", {}).get(name, {}))
            result = descriptor.build_plan(manifest=manifest, parameters=values, tasks=tuple(tasks))
            result = [result] if isinstance(result, TaskContract) else list(result)
            for task in result:
                if not isinstance(task, TaskContract) or task.capability != name:
                    raise ValueError("Extension plan factory returned an invalid capability contract")
                task.method = descriptor.method
            built.extend(result)
        return built

    def run_validator(
        self,
        contract: TaskContract,
        result: TaskResult,
        artifact_registry: ArtifactRegistry,
    ) -> Any:
        """Run a descriptor validator as an additional independent audit hook."""

        descriptor = self.describe(contract.capability, result.method_used or contract.method)
        return descriptor.run_validator(contract, result, artifact_registry)

    def extract_evidence(
        self,
        contract: TaskContract,
        result: TaskResult,
        report: Any,
        artifact_registry: ArtifactRegistry,
    ) -> list[Any]:
        """Run a descriptor evidence hook only after an accepted audit.

        Returned nodes are checked at this boundary so a plugin cannot attach
        another task's artifact or reuse an evidence ID.  The orchestrator
        remains responsible for adding accepted nodes to its EvidenceGraph.
        """

        if report is None:
            raise ValueError("Descriptor evidence extraction requires the independent audit report")
        if not report.overall_passed or report.stop_rule_triggered:
            return []
        descriptor = self.describe(contract.capability, result.method_used or contract.method)
        raw_nodes = descriptor.extract_evidence(contract, result, report, artifact_registry)
        if raw_nodes is None:
            return []
        if isinstance(raw_nodes, (str, bytes)):
            raise TypeError("evidence_extractor must return EvidenceNode objects, not text")
        try:
            nodes = list(raw_nodes)
        except TypeError as exc:
            nodes = [raw_nodes]
        seen_ids: set[str] = set()
        output_uris = set(result.output_artifacts)
        validated = []
        for node in nodes:
            if not hasattr(node, "evidence_id"):
                raise TypeError("evidence_extractor must return EvidenceNode objects")
            evidence_id = str(node.evidence_id)
            if not evidence_id or evidence_id in seen_ids:
                raise ValueError(f"Evidence extractor returned duplicate/empty evidence ID: {evidence_id!r}")
            if node.source_task_id != contract.task_id:
                raise ValueError("Descriptor evidence must cite the current task as source_task_id")
            if not node.audit_passed:
                raise ValueError("Descriptor evidence must be marked audit_passed after independent audit")
            if not node.source_artifact_uris or not set(node.source_artifact_uris).issubset(output_uris):
                raise ValueError("Descriptor evidence must cite only current task output artifacts")
            if any(not artifact_registry.exists(uri) for uri in node.source_artifact_uris):
                raise ValueError("Descriptor evidence cites an unregistered output artifact")
            seen_ids.add(evidence_id)
            validated.append(node)
        return validated

    def list_capabilities(self) -> Dict[str, List[Dict[str, Any]]]:
        result = {}
        for cap_name, impls in self._capabilities.items():
            result[cap_name] = [
                {
                    "implementation_id": cap.implementation_id,
                    "type": cap.implementation_type.value,
                    "suitable_for": cap.suitable_for,
                    "accepts_modalities": cap.accepts_modalities,
                    **{
                        "method": cap.get_descriptor().method,
                        "input_types": [value.value for value in cap.get_descriptor().input_types],
                        "output_types": [value.value for value in cap.get_descriptor().output_types],
                        "required_audit_ids": list(cap.get_descriptor().required_audit_ids),
                    },
                }
                for cap in impls.values()
            ]
        return result

    def execute_contract(self, contract: TaskContract, registry: ArtifactRegistry) -> TaskResult:
        """Executes a TaskContract through the registered capability with side-effect validation."""
        try:
            return self._execute_contract(contract, registry)
        except Exception as exc:
            return TaskResult(task_id=contract.task_id, capability=contract.capability,
                method_used=contract.method or "unresolved", status=TaskStatus.EXECUTION_FAILURE,
                error_type=ExecutionFailureType.DEPENDENCY_ERROR if isinstance(exc, ImportError) else ExecutionFailureType.CODE_ERROR,
                error_message=f"{type(exc).__name__}: {exc}")

    def _execute_contract(self, contract, registry):
        contract = self.validate_contract(contract, registry)
        capability = self.get(contract.capability, contract.method)

        # Load input payloads for side-effect pre-check
        input_payloads = {}
        for in_uri in contract.input_artifacts:
            input_payloads[in_uri] = registry.load_payload(in_uri)

        # Run capability execution
        task_result = capability.execute(contract, registry)

        # If execution succeeded, run side-effect validation against contract
        if task_result.status == TaskStatus.SUCCESS:
            missing = [u for u in contract.expected_outputs if u not in task_result.output_artifacts or not registry.exists(u)]
            if missing:
                task_result.status = TaskStatus.EXECUTION_FAILURE
                task_result.error_type = ExecutionFailureType.CODE_ERROR
                task_result.error_message = f"Required output artifacts missing: {missing}"
                return task_result
            output_payloads = {}
            for out_uri in task_result.output_artifacts:
                if registry.exists(out_uri):
                    descriptor = capability.get_descriptor()
                    if descriptor.validate_types and descriptor.output_types and registry.get_metadata(out_uri).type not in descriptor.output_types:
                        raise ValueError("Output artifact type violates capability descriptor")
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
