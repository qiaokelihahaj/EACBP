"""Capability interfaces and the shared capability declaration.

Historically the execution classes carried most of their contract metadata as
ad-hoc attributes.  That made it very easy for a planner, registry, and
auditor to disagree about a method.  :class:`CapabilityDescriptor` is the
small, executable declaration used at those boundaries now.  It deliberately
does not replace ``BaseCapability``: an implementation still owns execution,
while the descriptor owns the method contract and optional integration hooks.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import asdict, is_dataclass
from typing import Callable, Iterable, List, Dict, Any, Optional, Sequence, Tuple
from enum import Enum

from pydantic import BaseModel
from eacbp.schemas.task import TaskContract, TaskResult
from eacbp.schemas.artifact import ArtifactType
from eacbp.artifact.registry import ArtifactRegistry


class ImplementationType(str, Enum):
    PYTHON_TOOL = "python_tool"
    R_TOOL = "r_tool"
    AGENT_ADAPTER = "agent_adapter"
    CONTAINER = "container"


def _normalise_types(values: Optional[Sequence[ArtifactType]]) -> Tuple[ArtifactType, ...]:
    """Return artifact types in a stable, duplicate-free tuple.

    Descriptor values are often declared with a list for readability.  A
    tuple makes the declaration safe to share between capability instances,
    and preserving enum values keeps comparisons against artifact metadata
    unambiguous.
    """

    result: list[ArtifactType] = []
    for value in values or ():
        if not isinstance(value, ArtifactType):
            try:
                value = ArtifactType(value)
            except (TypeError, ValueError) as exc:
                raise ValueError(f"Unknown artifact type in capability descriptor: {value!r}") from exc
        if value not in result:
            result.append(value)
    return tuple(result)


def _invoke_hook(hook: Callable[..., Any], *args: Any) -> Any:
    """Invoke a descriptor hook using its documented positional protocol.

    Validators receive ``(contract, result, registry)`` and evidence
    extractors receive ``(contract, result, report, registry)``.  The
    descriptor boundary deliberately does not guess a callable's arity or
    rewrite keyword-only parameters: doing so can hide a real plugin bug or a
    ``TypeError`` raised by the hook itself.  ``audit`` remains a small
    compatibility escape hatch for validator objects used by older plugins;
    the selected method still receives the complete canonical argument list.
    """

    target = getattr(hook, "audit", hook)
    if not callable(target):
        raise TypeError("descriptor hook must be callable")
    return target(*args)


class CapabilityDescriptor:
    """Executable declaration shared by planning, execution, and auditing.

    ``validator`` and ``evidence_extractor`` are integration hooks, not a
    replacement for the independent scientific auditor.  The orchestrator
    invokes a validator as an additional audit and only invokes an evidence
    extractor after that audit has passed.  ``plan_factory`` is used by
    extension planners and is optional for ordinary capabilities.

    The public constructor intentionally mirrors the existing vocabulary in
    ``BaseCapability`` (``capability_name`` and ``method``) so declarations
    remain explicit and discoverable.
    """

    __slots__ = (
        "capability_name", "method", "parameter_model", "input_types", "output_types",
        "required_audit_ids", "validator", "evidence_extractor", "plan_factory",
        "method_aliases", "description", "validate_types",
    )

    def __init__(
        self,
        capability_name: str,
        method: str,
        parameter_model: Any = None,
        input_types: Optional[Sequence[ArtifactType]] = None,
        output_types: Optional[Sequence[ArtifactType]] = None,
        required_audit_ids: Optional[Sequence[str]] = None,
        validator: Optional[Callable[..., Any]] = None,
        evidence_extractor: Optional[Callable[..., Any]] = None,
        plan_factory: Optional[Callable[..., Any]] = None,
        method_aliases: Optional[Sequence[str]] = None,
        description: str = "",
        validate_types: bool = True,
    ) -> None:
        if not isinstance(capability_name, str) or not capability_name.strip():
            raise ValueError("Capability descriptor requires a non-empty capability_name")
        if not isinstance(method, str) or not method.strip():
            raise ValueError("Capability descriptor requires a non-empty method")
        self.capability_name = capability_name.strip()
        self.method = method.strip()
        self.parameter_model = parameter_model
        self.input_types = _normalise_types(input_types)
        self.output_types = _normalise_types(output_types)
        self.required_audit_ids = tuple(dict.fromkeys(str(value) for value in (required_audit_ids or ()) if str(value).strip()))
        self.validator = validator
        self.evidence_extractor = evidence_extractor
        self.plan_factory = plan_factory
        self.method_aliases = tuple(dict.fromkeys(str(value) for value in (method_aliases or ()) if str(value).strip()))
        self.description = str(description or "")
        self.validate_types = bool(validate_types)

    @property
    def capability(self) -> str:
        """Short alias useful when serialising a descriptor for discovery."""

        return self.capability_name

    @property
    def required_audits(self) -> Tuple[str, ...]:
        """Compatibility/readability alias for ``required_audit_ids``."""

        return self.required_audit_ids

    @classmethod
    def from_capability(cls, capability: "BaseCapability") -> "CapabilityDescriptor":
        """Build a declaration for legacy capabilities with no explicit one."""

        return cls(
            capability_name=capability.capability_name,
            method=capability.implementation_id,
            input_types=capability.accepts_types,
            output_types=capability.output_types,
            required_audit_ids=getattr(capability, "required_audit_ids", None)
            or getattr(capability, "required_audit_checks", None),
            method_aliases=tuple(getattr(capability, "legacy_aliases", {}).keys()),
            validate_types=False,
        )

    def validate_parameters(self, parameters: Optional[Dict[str, Any]]) -> Dict[str, Any]:
        """Validate and normalise parameters while retaining orchestration extras.

        A descriptor may use a Pydantic model, a dataclass, or a callable
        returning a mapping/model.  Parameters not represented by a Pydantic
        model are merged back into the result so metadata injected by the
        orchestrator (study id, target branch, FDR family, and similar
        routing context) is never silently discarded.
        """

        if parameters is None:
            raw: Dict[str, Any] = {}
        elif isinstance(parameters, dict):
            raw = dict(parameters)
        else:
            try:
                raw = dict(parameters)
            except (TypeError, ValueError) as exc:
                raise TypeError("Capability parameters must be a mapping") from exc
        model = self.parameter_model
        if model is None:
            return raw

        if isinstance(model, type) and issubclass(model, BaseModel):
            validated = model.model_validate(raw)
            dumped = validated.model_dump(mode="python")
        elif hasattr(model, "model_validate"):
            validated = model.model_validate(raw)
            dumped = validated.model_dump(mode="python") if hasattr(validated, "model_dump") else dict(validated)
        elif isinstance(model, type) and is_dataclass(model):
            validated = model(**raw)
            dumped = asdict(validated)
        elif callable(model):
            validated = model(raw)
            if isinstance(validated, BaseModel):
                dumped = validated.model_dump(mode="python")
            elif is_dataclass(validated):
                dumped = asdict(validated)
            elif isinstance(validated, dict):
                dumped = dict(validated)
            else:
                raise TypeError("Capability parameter model must return a mapping or model")
        else:
            raise TypeError("parameter_model must be a Pydantic model, dataclass, or callable")
        if not isinstance(dumped, dict):
            raise TypeError("Capability parameter model did not produce a mapping")
        # Preserve extension/routing values that a permissive model did not
        # declare.  Models configured with extra='forbid' still reject them.
        return {**raw, **dumped}

    def run_validator(self, contract: Any, result: Any, registry: Any) -> Any:
        """Execute the optional independent validator hook, if declared."""

        if self.validator is None:
            return None
        return _invoke_hook(self.validator, contract, result, registry)

    def extract_evidence(self, contract: Any, result: Any, report: Any, registry: Any) -> Any:
        """Execute the optional evidence hook for an already audited result."""

        if self.evidence_extractor is None:
            return []
        return _invoke_hook(self.evidence_extractor, contract, result, report, registry)

    def build_plan(self, *args: Any, **kwargs: Any) -> Any:
        """Run the optional plan factory used by extension planners."""

        if self.plan_factory is None:
            return None
        return self.plan_factory(*args, **kwargs)

    def to_dict(self) -> Dict[str, Any]:
        """Return a JSON-friendly discovery summary (hooks are intentionally opaque)."""

        return {
            "capability_name": self.capability_name,
            "method": self.method,
            "parameter_model": getattr(self.parameter_model, "__name__", None)
            or (self.parameter_model.__class__.__name__ if self.parameter_model is not None else None),
            "input_types": [value.value for value in self.input_types],
            "output_types": [value.value for value in self.output_types],
            "required_audit_ids": list(self.required_audit_ids),
            "method_aliases": list(self.method_aliases),
            "has_validator": self.validator is not None,
            "has_evidence_extractor": self.evidence_extractor is not None,
            "has_plan_factory": self.plan_factory is not None,
            "validate_types": self.validate_types,
            "description": self.description,
        }


class BaseCapability(ABC):
    """Abstract interface for all computational capabilities (tools, pipelines, agent adapters)."""

    def __init__(
        self,
        capability_name: str,
        implementation_id: str,
        implementation_type: ImplementationType = ImplementationType.PYTHON_TOOL,
        accepts_modalities: Optional[List[str]] = None,
        accepts_types: Optional[List[ArtifactType]] = None,
        requires_keys: Optional[List[str]] = None,
        suitable_for: Optional[List[str]] = None,
        output_types: Optional[List[ArtifactType]] = None,
        descriptor: Optional[CapabilityDescriptor] = None,
    ):
        self.capability_name = capability_name
        self.implementation_id = implementation_id
        self.implementation_type = implementation_type
        self.accepts_modalities = accepts_modalities or ["scRNA"]
        self.accepts_types = accepts_types or [ArtifactType.ANNDATA]
        self.requires_keys = requires_keys or []
        self.suitable_for = suitable_for or []
        self.output_types = output_types or [ArtifactType.ANNDATA]
        # Legacy capabilities get a descriptor lazily from their existing
        # metadata. Explicit declarations (including plugin capabilities) set
        # this attribute in their constructor.
        self.descriptor = descriptor

    def get_descriptor(self) -> CapabilityDescriptor:
        """Return the shared declaration for this implementation.

        Advanced extension descriptors live in a separate module to keep the
        base class importable without importing optional scientific packages.
        A legacy capability transparently falls back to its constructor
        metadata, so this change is backwards compatible.
        """

        if self.descriptor is not None:
            return self.descriptor
        try:
            from eacbp.capabilities.advanced_descriptors import descriptor_for_capability
            descriptor = descriptor_for_capability(self.capability_name, self.implementation_id)
        except ImportError:
            descriptor = None
        if descriptor is not None:
            self.descriptor = descriptor
            return descriptor
        descriptor = CapabilityDescriptor.from_capability(self)
        self.descriptor = descriptor
        return descriptor

    @abstractmethod
    def execute(self, contract: TaskContract, registry: ArtifactRegistry) -> TaskResult:
        """Executes the contracted computation, registering new artifacts in registry and returning TaskResult."""
        pass
