"""Typed, explicit provenance and design contract for a research acceptance run."""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Literal, Optional

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator


class StrictModel(BaseModel):
    model_config = ConfigDict(extra="forbid")


class ContrastManifest(StrictModel):
    column: str = Field(min_length=1)
    tested: str = Field(min_length=1, description="Numerator / effect direction")
    reference: str = Field(min_length=1, description="Denominator / reference level")

    @model_validator(mode="after")
    def different_levels(self):
        if self.tested == self.reference:
            raise ValueError("contrast tested and reference levels must differ")
        return self


class DonorManifest(StrictModel):
    column: str = Field(min_length=1)
    biological_unit: str = Field(min_length=1)
    paired: bool
    minimum_per_condition: int = Field(default=2, ge=2)


class ReferenceEnvironment(StrictModel):
    python_version: str = Field(min_length=1)
    packages: dict[str, str] = Field(min_length=1)
    lockfile_path: Optional[str] = None
    lockfile_sha256: Optional[str] = None

    @field_validator("lockfile_sha256")
    @classmethod
    def valid_lock_hash(cls, value):
        if value is not None and not re.fullmatch(r"[0-9a-fA-F]{64}", value):
            raise ValueError("lockfile_sha256 must be a 64-character SHA256 hex digest")
        return value.lower() if value else value

    @model_validator(mode="after")
    def python_version_has_one_source(self):
        reserved = {"python", "python_version"}
        if any(name.strip().lower().replace("-", "_") in reserved for name in self.packages):
            raise ValueError("declare Python only through python_version, not packages")
        return self


class ResourceLimits(StrictModel):
    max_cells: Optional[int] = Field(default=None, ge=1)
    max_genes: Optional[int] = Field(default=None, ge=1)
    max_nonzero_counts: Optional[int] = Field(default=None, ge=1)
    max_input_bytes: Optional[int] = Field(default=None, ge=1)
    max_pseudobulk_working_bytes: Optional[int] = Field(default=None, ge=1)


class StatisticsSettings(StrictModel):
    alpha: float = Field(default=0.05, gt=0, lt=1, allow_inf_nan=False)
    cooks_filter: bool = True
    independent_filter: bool = True
    fit_type: Literal["parametric", "mean"] = "parametric"
    size_factors_fit_type: Literal["ratio", "poscounts", "iterative"] = "ratio"
    rtol: float = Field(default=1e-6, ge=0, allow_inf_nan=False)
    atol: float = Field(default=1e-8, ge=0, allow_inf_nan=False)


class AcceptanceManifest(StrictModel):
    schema_version: Literal[1] = 1
    study_id: str = Field(min_length=1, pattern=r"^[A-Za-z0-9_-]+$")
    dataset_kind: Literal["synthetic", "real"]
    data_path: str = Field(min_length=1)
    sha256: str = Field(pattern=r"^[0-9a-fA-F]{64}$")
    source: str = Field(min_length=1, description="Accession or complete citation")
    source_url: Optional[str] = None
    license: str = Field(min_length=1)
    species: str = Field(min_length=1)
    species_metadata_key: Optional[str] = None
    counts_layer: str = Field(min_length=1)
    contrast: ContrastManifest
    donor: DonorManifest
    batch_column: Optional[str] = None
    covariates: list[str] = Field(default_factory=list)
    design_formula: Optional[str] = None
    statistics: StatisticsSettings = Field(default_factory=StatisticsSettings)
    resource_limits: ResourceLimits = Field(default_factory=ResourceLimits)
    reference_environment: Optional[ReferenceEnvironment] = None

    @field_validator("sha256")
    @classmethod
    def normalize_sha256(cls, value):
        return value.lower()

    @model_validator(mode="after")
    def validate_manifest_relationships(self):
        if self.batch_column and self.batch_column not in self.covariates:
            self.covariates.append(self.batch_column)
        if len(set(self.covariates)) != len(self.covariates):
            raise ValueError("covariate names must be unique")
        if self.contrast.column == self.donor.column:
            raise ValueError("condition and donor metadata columns must differ")
        if self.dataset_kind == "real":
            if self.reference_environment is None:
                raise ValueError("real research acceptance requires a fixed reference_environment")
            if not _filled(self.source) or "replace" in self.source.strip().lower():
                raise ValueError("real research acceptance requires a source accession or complete citation")
            if not self.source_url or not _filled(self.source_url) or "replace" in self.source_url.strip().lower():
                raise ValueError("real research acceptance requires a source_url or accession link")
            if not _filled(self.license) or "replace" in self.license.strip().lower():
                raise ValueError("real research acceptance requires a known license")
            if not _filled(self.species) or "replace" in self.species.strip().lower():
                raise ValueError("real research acceptance requires a declared species")
            if not self.reference_environment.lockfile_path or not self.reference_environment.lockfile_sha256:
                raise ValueError("real research acceptance requires a hashed environment lockfile")
        return self


def _filled(value: Optional[str]) -> bool:
    if value is None:
        return False
    return value.strip().lower() not in {"", "tbd", "unknown", "replace_me", "required"}


def load_manifest(path: str | Path) -> tuple[AcceptanceManifest, Path]:
    """Load a manifest and return it with its resolved parent directory."""
    path = Path(path).expanduser().resolve()
    payload = json.loads(path.read_text(encoding="utf-8"))
    return AcceptanceManifest.model_validate(payload), path.parent
