"""
Study Manifest schema defining biological, experimental, data, policy, and reproducibility specifications.
"""

from typing import List, Dict, Optional
from pydantic import BaseModel, Field


class BiologicalDesign(BaseModel):
    species: str = Field(..., description="Species name, e.g., mus_musculus, homo_sapiens")
    tissue: str = Field(..., description="Tissue or anatomical region, e.g., brain, cortex")
    disease: Optional[str] = Field(None, description="Disease context, e.g., Alzheimer, healthy")
    conditions: List[str] = Field(default_factory=list, description="Experimental conditions, e.g., ['AD', 'control']")
    target_cell_types: List[str] = Field(default_factory=list, description="Target cell types of interest, e.g., ['microglia']")


class ExperimentalDesign(BaseModel):
    biological_unit: str = Field("mouse", description="Biological replication unit, e.g. mouse, patient, donor")
    batches: List[str] = Field(default_factory=list, description="Batch identifiers")
    total_samples: int = Field(1, description="Total number of biological samples")
    donor_replicates_per_condition: Dict[str, int] = Field(
        default_factory=dict,
        description="Number of biological replicates per condition, e.g. {'AD': 6, 'control': 6}"
    )


class DataSpec(BaseModel):
    modalities: List[str] = Field(default_factory=lambda: ["scRNA"], description="Data modalities, e.g., ['scRNA', 'spatial', 'FASTQ']")
    raw_artifact_uri: Optional[str] = Field(None, description="URI of the initial raw data artifact")
    has_raw_fastq: bool = Field(False, description="Whether raw paired-end FASTQ reads are provided for quantification")
    fastq_paths: Dict[str, Dict[str, str]] = Field(default_factory=dict, description="Mapping of sample ID to R1/R2 fastq file paths")
    chemistry: str = Field("10xv3", description="Single-cell chemistry, e.g., 10xv3, 10xv2")
    has_spatial_coordinates: bool = Field(False, description="Whether spatial 2D/3D coordinates are present")
    has_rna_velocity: bool = Field(False, description="Whether spliced/unspliced velocity layers are present")
    has_multiple_timepoints: bool = Field(False, description="Whether temporal longitudinal series are present")


class Hypotheses(BaseModel):
    user_provided: List[str] = Field(default_factory=list, description="User-specified prior hypotheses")
    prior_generated: List[str] = Field(default_factory=list, description="Literature or model generated priors")


class Constraints(BaseModel):
    max_runtime_hours: float = Field(12.0, description="Max runtime constraint in hours")
    gpu_allowed: bool = Field(True, description="Whether GPU hardware acceleration is permitted")
    min_biological_replicates: int = Field(2, description="Minimum donor replicates required for confirmatory DEG")


class AnalysisPolicy(BaseModel):
    discovery_mode: bool = Field(True, description="Discovery mode does not bias initial feature selection with prior hypotheses")
    prior_guided_analysis: bool = Field(False, description="Whether analysis explicitly tests a user-guided prior axis")
    strict_reproducibility: bool = Field(True, description="Strict random seed & dependency provenance tracking")
    prefer_pseudobulk: bool = Field(True, description="Prefer pseudobulk differential expression when biological replicates >= 3")


class ReproducibilityConfig(BaseModel):
    random_seed: int = Field(42, description="Fixed random seed for deterministic pipelines")
    save_intermediate_artifacts: bool = Field(True, description="Persist all intermediate AnnData checkpoints")


class StudyManifest(BaseModel):
    study_id: str = Field(..., description="Unique study identifier, e.g., AD_mouse_001")
    title: str = Field("Single-cell Study", description="Human-readable title of the study")
    biological_design: BiologicalDesign
    experimental_design: ExperimentalDesign = Field(default_factory=ExperimentalDesign)
    data: DataSpec = Field(default_factory=DataSpec)
    hypotheses: Hypotheses = Field(default_factory=Hypotheses)
    constraints: Constraints = Field(default_factory=Constraints)
    analysis_policy: AnalysisPolicy = Field(default_factory=AnalysisPolicy)
    reproducibility: ReproducibilityConfig = Field(default_factory=ReproducibilityConfig)
