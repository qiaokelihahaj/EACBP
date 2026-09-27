"""Static parameter contracts for the highest-cost/common built-in methods.

These models validate supplied options without injecting defaults into legacy
task signatures. Data-dependent eligibility still belongs to execution and the
independent scientific auditor. Routing metadata is an explicit shared base;
misspelled algorithm options are not silently accepted as metadata.
"""

from pathlib import Path
from typing import Annotated, Any, ClassVar, Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator


NonEmptyString = Annotated[str, Field(min_length=1)]
NonNegativeInt = Annotated[int, Field(ge=0)]
PositiveInt = Annotated[int, Field(ge=1)]
ReplicateCount = Annotated[int, Field(ge=2)]
PositiveFloat = Annotated[float, Field(gt=0)]
Probability = Annotated[float, Field(gt=0, lt=1)]


class MethodParameters(BaseModel):
    model_config = ConfigDict(extra="forbid", strict=True, allow_inf_nan=False)
    preserve_omitted_parameters: ClassVar[bool] = True

    study_id: NonEmptyString | None = None
    target_cell_type: NonEmptyString | None = None
    target_branch: NonEmptyString | None = None
    target_dependent: bool | None = None
    target_provenance: dict[str, Any] | None = None
    fdr_family: NonEmptyString | None = None
    species: NonEmptyString | None = None
    tissue: NonEmptyString | None = None
    random_seed: NonNegativeInt = 42
    external_resource_sha256: dict[str, str] | None = None


class QCParameters(MethodParameters):
    min_genes: NonNegativeInt = 100
    max_mito_pct: Annotated[float, Field(ge=0, le=100)] = 20.0


class NormalizationParameters(MethodParameters):
    target_sum: PositiveFloat = 10000.0
    n_top_genes: NonNegativeInt = 300
    input_layer: NonEmptyString = "counts"


class ContrastParameters(MethodParameters):
    condition_col: NonEmptyString = "condition"
    donor_col: NonEmptyString | None = None
    # The existing form carries this metadata label with its contrast. It
    # does not implicitly add a batch covariate to a statistical design.
    batch_col: NonEmptyString | None = None
    condition_a: NonEmptyString | None = None
    condition_b: NonEmptyString | None = None
    min_replicates: ReplicateCount = 2
    allow_x_as_counts: bool = False

    @model_validator(mode="after")
    def validate_contrast(self):
        if (self.condition_a is None) != (self.condition_b is None):
            raise ValueError("condition_a and condition_b must be supplied together")
        if self.condition_a is not None and self.condition_a == self.condition_b:
            raise ValueError("condition_a and condition_b must be different")
        return self


class WelchDEGParameters(ContrastParameters):
    require_pseudobulk: bool = False
    pseudocount: PositiveFloat = 0.5


class PyDESeq2Parameters(ContrastParameters):
    min_donors: ReplicateCount = 2
    paired: bool = False
    paired_design: bool = False
    covariates: str | list[str] | tuple[str, ...] | None = None
    covariate_cols: str | list[str] | tuple[str, ...] | None = None
    counts_layer: NonEmptyString = "counts"
    raw_counts_layer: NonEmptyString = "counts"
    design_formula: NonEmptyString | None = None
    design: NonEmptyString | None = None
    n_cpus: PositiveInt | None = 1
    n_processes: PositiveInt | None = 1
    low_memory: bool = False
    fit_type: Literal["parametric", "mean"] = "parametric"
    size_factors_fit_type: Literal["ratio", "poscounts", "iterative"] = "ratio"
    alpha: Probability = 0.05
    cooks_filter: bool = True
    independent_filter: bool = True
    contrast: list[str] | list[float] | tuple[str, ...] | tuple[float, ...] | None = None
    contrast_label: NonEmptyString | None = None
    contrast_effect_definition: NonEmptyString | None = None

    @field_validator("contrast", mode="before")
    @classmethod
    def normalize_array_contrast(cls, value):
        # The Python API also accepts a one-dimensional NumPy contrast array.
        return value.tolist() if hasattr(value, "tolist") else value


class CellBenderParameters(MethodParameters):
    unfiltered_input_path: NonEmptyString | Path
    executable: NonEmptyString | Path
    output_path: NonEmptyString | Path
    run_cwd: NonEmptyString | Path | None = None
    subcommand: Literal["remove-background"] = "remove-background"
    allow_overwrite: bool = False
    extra_args: list[str] | tuple[str, ...] | None = None
    cli_args: list[str] | tuple[str, ...] | None = None
    timeout_sec: PositiveFloat = 3600.0
    cell_id_key: NonEmptyString = "cell_id"
    gene_id_key: NonEmptyString = "gene_name"
    report_path: NonEmptyString | Path | None = None
    metrics_path: NonEmptyString | Path | None = None
    log_path: NonEmptyString | Path | None = None

    @field_validator("extra_args", "cli_args")
    @classmethod
    def protect_io_arguments(cls, value):
        if any(arg in {"--input", "--output"} or arg.startswith(("--input=", "--output="))
               for arg in value or ()):
            raise ValueError("extra_args/cli_args cannot override explicit --input or --output")
        return value

    @field_validator("output_path")
    @classmethod
    def require_h5_output(cls, value):
        if Path(value).suffix.casefold() != ".h5":
            raise ValueError("CellBender output_path must end with .h5")
        return value
