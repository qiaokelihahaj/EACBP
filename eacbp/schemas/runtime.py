"""Typed runtime state used by the scientific orchestrator.

The public orchestrator API historically accepted a ``dict`` called
``current_state``.  The dictionary mixed three different kinds of values:
user supplied run options, observations learned while executing the DAG, and
bookkeeping controls such as ``resume``.  This module keeps that API as a
compatibility boundary while giving the execution loop explicit containers
for each concern.

``RunConfig`` owns values supplied by the caller.  ``PlanningContext`` owns
derived and observed values.  ``ExecutionState`` is the mutable coordinator
that exposes a legacy mapping to the planner/router.  Configuration wins when
the same key exists in both places, so a result metric can no longer silently
change a user-selected method or extension.
"""

from __future__ import annotations

from copy import deepcopy
from typing import Any, ClassVar, Dict, Iterable, List, Literal, Mapping, Optional, Set, Union

from pydantic import BaseModel, ConfigDict, Field, StrictBool


class RunConfig(BaseModel):
    """Caller supplied analysis options.

    The fields below cover the options used by the built-in planner.  Extra
    values remain accepted deliberately: external capability plugins have
    historically added their own options to ``current_state``.  Unknown
    options are retained as config values and therefore receive the same
    precedence guarantee as declared fields.
    """

    model_config = ConfigDict(extra="allow", validate_assignment=True)

    method_profile: Optional[Literal["baseline", "standard"]] = None
    mode: Optional[Literal["demo", "real"]] = None
    resume: Optional[StrictBool] = None
    method_overrides: Optional[Dict[str, str]] = None
    capability_parameters: Optional[Dict[str, Dict[str, Any]]] = None
    target_parameters: Optional[Dict[str, Dict[str, Dict[str, Any]]]] = None
    analysis_extensions: Optional[Dict[str, Union[StrictBool, Dict[str, Any]]]] = None
    advanced_analysis: Optional[StrictBool] = None

    quant_tool: Optional[str] = None
    target_gene: Optional[str] = None
    target_genes: Optional[List[str]] = None
    perturbation_type: Optional[str] = None
    perturbation_targets: Optional[List[str]] = None
    cellrank_terminal_states: Optional[Dict[str, List[str]]] = None

    threads: Optional[int] = None
    n_cells: Optional[int] = None
    n_genes: Optional[int] = None
    work_dir: Optional[str] = None
    index_path: Optional[str] = None
    t2g_path: Optional[str] = None
    star_bin: Optional[str] = None
    genome_dir: Optional[str] = None
    whitelist_path: Optional[str] = None
    gtf_path: Optional[str] = None
    num_reads: Optional[int] = None

    has_fastq: Optional[StrictBool] = None
    include_spatial: Optional[StrictBool] = None
    include_cci: Optional[StrictBool] = None
    run_cell_communication: Optional[StrictBool] = None
    include_adapters: Optional[StrictBool] = None
    adapters: Optional[List[str]] = None
    run_spacell_adapter: Optional[StrictBool] = None
    run_geneagent_adapter: Optional[StrictBool] = None
    run_chatcell_adapter: Optional[StrictBool] = None
    include_knowledge: Optional[StrictBool] = None
    run_knowledge_engine: Optional[StrictBool] = None
    full_e2e: Optional[StrictBool] = None
    run_perturbation: Optional[StrictBool] = None
    include_perturbation: Optional[StrictBool] = None
    run_compound_perturbation: Optional[StrictBool] = None
    prior_guided: Optional[StrictBool] = None
    batch_effect_possible: Optional[StrictBool] = None

    # A metric with one of these names must never become a new routing or
    # execution option merely because a capability returned it.  If the caller
    # supplied the key, its actual value is still restored by ``as_mapping``.
    # Observation-only keys such as n_cells/n_genes are intentionally omitted:
    # dataset audit results use those to adapt later tasks.
    PROTECTED_KEYS: ClassVar[Set[str]] = {
        "method_profile", "mode", "resume", "method_overrides",
        "capability_parameters", "target_parameters", "analysis_extensions", "advanced_analysis",
        "quant_tool", "target_gene", "target_genes", "perturbation_type",
        "perturbation_targets", "cellrank_terminal_states", "threads",
        "work_dir", "index_path", "t2g_path", "star_bin", "genome_dir",
        "whitelist_path", "gtf_path", "num_reads", "has_fastq",
        "include_spatial", "include_cci", "run_cell_communication",
        "include_adapters", "adapters", "run_spacell_adapter",
        "run_geneagent_adapter", "run_chatcell_adapter", "include_knowledge",
        "run_knowledge_engine", "full_e2e", "run_perturbation",
        "include_perturbation", "run_compound_perturbation", "prior_guided",
    }

    def to_mapping(self) -> Dict[str, Any]:
        """Return a deep copied mapping without unset optional defaults."""

        return deepcopy(self.model_dump(exclude_none=True, exclude_unset=True))

    @classmethod
    def from_mapping(cls, value: Optional[Mapping[str, Any] | "RunConfig"]) -> "RunConfig":
        """Coerce the legacy mapping or an existing config into a copy."""

        if value is None:
            return cls()
        if isinstance(value, cls):
            return value.model_copy(deep=True)
        if not isinstance(value, Mapping):
            raise TypeError("run configuration must be a mapping or RunConfig")
        return cls.model_validate(deepcopy(dict(value)))

    def get(self, key: str, default: Any = None) -> Any:
        """Mapping-style access retained for router/planner compatibility."""

        return self.to_mapping().get(key, default)

    def __contains__(self, key: str) -> bool:
        return key in self.to_mapping()

    def __getitem__(self, key: str) -> Any:
        return self.to_mapping()[key]


class PlanningContext(BaseModel):
    """Values derived from the manifest and observations made by the DAG."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    derived: Dict[str, Any] = Field(default_factory=dict)
    observed: Dict[str, Any] = Field(default_factory=dict)
    branches: Dict[str, Dict[str, Any]] = Field(default_factory=dict)
    decisions: List[Dict[str, Any]] = Field(default_factory=list)

    # Dataset/QC audit outputs are observations rather than user-selected
    # execution options.  Keeping their legacy seed values here allows a
    # later audit to replace stale metadata such as ``batches``.
    OBSERVATION_KEYS: ClassVar[Set[str]] = {
        "batches",
        "observed_conditions",
        "condition_metadata_complete",
        "min_replicates",
        "biological_replication_sufficient",
        "paired_donors_observed",
    }

    def set_derived(self, key: str, value: Any) -> None:
        self.derived[key] = deepcopy(value)

    def update_observed(self, values: Mapping[str, Any]) -> None:
        """Merge metrics into the observation plane without touching config."""

        self.observed.update(deepcopy(dict(values)))

    def add_decisions(self, values: Iterable[Mapping[str, Any]]) -> None:
        self.decisions.extend(deepcopy([dict(item) for item in values]))

    def as_mapping(self, config: RunConfig) -> Dict[str, Any]:
        """Build the state view expected by legacy planner/router code.

        Observations are intentionally applied first.  Derived values then
        override an observation of the same name, and explicit caller config
        has final precedence.
        """

        result = {
            key: deepcopy(value)
            for key, value in self.observed.items()
            if key not in config.PROTECTED_KEYS
        }
        # Derived values are produced by the planner itself (for example the
        # inferred method profile), so they are trusted and remain visible.
        # Explicit config below still wins over a derived default.
        result.update(deepcopy(self.derived))
        result.update(config.to_mapping())
        return result


class ExecutionState(BaseModel):
    """Mutable, typed coordinator state for one study execution."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    config: RunConfig = Field(default_factory=RunConfig)
    planning: PlanningContext = Field(default_factory=PlanningContext)
    resume_requested: bool = False
    completed: Dict[str, str] = Field(default_factory=dict)

    @classmethod
    def from_input(cls, value: Optional[Mapping[str, Any] | RunConfig]) -> "ExecutionState":
        # Validate controls before splitting them out of the legacy mapping.
        # In particular, bool("false") must never request a resumed run.
        raw = RunConfig.from_mapping(value).to_mapping()
        resume = bool(raw.pop("resume", False))
        observed = {
            key: raw.pop(key)
            for key in PlanningContext.OBSERVATION_KEYS
            if key in raw
        }
        return cls(
            config=RunConfig.from_mapping(raw),
            planning=PlanningContext(observed=observed),
            resume_requested=resume,
        )

    @property
    def current_state(self) -> Dict[str, Any]:
        """Legacy planner/router view; always returns an independent copy."""

        return self.planning.as_mapping(self.config)

    def set_derived(self, key: str, value: Any) -> None:
        self.planning.set_derived(key, value)

    def record_metrics(self, metrics: Mapping[str, Any], branch: Optional[str] = None) -> None:
        if branch is None:
            self.planning.update_observed(metrics)
        else:
            self.planning.branches.setdefault(branch, {}).update(deepcopy(dict(metrics)))

    def state_for_branch(self, branch: Optional[str]) -> Dict[str, Any]:
        """Expose shared observations plus only this target's observations."""
        return self.preview_with_metrics(self.planning.branches.get(branch, {})) if branch else self.current_state

    def preview_with_metrics(self, metrics: Mapping[str, Any]) -> Dict[str, Any]:
        """Return the next planner view without mutating this state."""

        context = self.planning.model_copy(deep=True)
        context.update_observed(metrics)
        return context.as_mapping(self.config)

    def record_decisions(self, decisions: Iterable[Mapping[str, Any]]) -> None:
        self.planning.add_decisions(decisions)

    def mark_completed(self, task_id: str, status: Any) -> None:
        value = getattr(status, "value", status)
        self.completed[task_id] = str(value)

    def config_keys(self) -> Set[str]:
        return set(self.config.to_mapping())

    def to_mapping(self) -> Dict[str, Any]:
        return self.current_state
