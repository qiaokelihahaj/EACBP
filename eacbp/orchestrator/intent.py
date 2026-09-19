"""
Intent parsing for a study manifest.

Parsing is intentionally conservative: a short substring such as ``ad`` or
``ko`` is not evidence of a disease or intervention, and sample counts are
only populated when an explicit sample manifest is supplied.
"""

from __future__ import annotations

import re
from collections import Counter
from typing import Any, Dict, List, Mapping, Optional

from eacbp.schemas.study import (
    AnalysisPolicy,
    BiologicalDesign,
    DataSpec,
    ExperimentalDesign,
    Hypotheses,
    ReproducibilityConfig,
    StudyManifest,
)


def _has_token(text: str, token: str) -> bool:
    """Match a standalone Latin token, avoiding substring false positives."""

    return bool(re.search(rf"(?<![A-Za-z0-9]){re.escape(token)}(?![A-Za-z0-9])", text, re.IGNORECASE))


def _has_any_token(text: str, tokens: List[str]) -> bool:
    return any(_has_token(text, token) for token in tokens)


def _has_any_phrase(text: str, phrases: List[str]) -> bool:
    return any(phrase.lower() in text for phrase in phrases)


def _has_ko(text: str) -> bool:
    return (
        _has_any_token(text, ["cko", "ko"])
        or _has_any_phrase(
            text,
            [
                "knockout",
                "knock-out",
                "conditional knockout",
                "条件性敲除",
                "敲除",
            ],
        )
    )


def _is_kat8(text: str) -> bool:
    return _has_any_token(text, ["kat8", "mof", "myst1"]) or _has_any_phrase(
        text,
        ["h4k16ac", "组蛋白乙酰化"],
    )


def _is_ad(text: str) -> bool:
    # ``ad`` must be a complete abbreviation.  In particular, it does not
    # match the beginning of ``adrenal``.
    return _has_any_token(text, ["ad", "alzheimer"]) or _has_any_phrase(
        text,
        ["阿尔茨海默", "痴呆"],
    )


def _extract_sample_manifest(sample_manifest: Optional[Mapping[str, Any]]) -> Dict[str, Any]:
    """Extract explicit sample metadata without inventing missing values."""

    if not isinstance(sample_manifest, Mapping):
        return {
            "paths": {},
            "total_samples": 0,
            "replicates": {},
            "batches": [],
            "conditions": [],
        }
    raw_samples = sample_manifest.get("samples", sample_manifest)
    if not isinstance(raw_samples, Mapping):
        return {
            "paths": {},
            "total_samples": 0,
            "replicates": {},
            "batches": [],
            "conditions": [],
        }

    paths: Dict[str, Dict[str, str]] = {}
    conditions: List[str] = []
    batches: List[str] = []
    replicate_counts: Counter[str] = Counter()
    for name, spec in raw_samples.items():
        if not isinstance(spec, Mapping):
            continue
        metadata = spec.get("metadata", {})
        metadata = metadata if isinstance(metadata, Mapping) else {}
        merged = dict(spec)
        merged.update(metadata)
        condition = str(merged.get("condition", "")).strip()
        donor = str(merged.get("donor", "")).strip()
        batch = str(merged.get("batch", "")).strip()
        if condition:
            conditions.append(condition)
        if batch:
            batches.append(batch)
        if condition and donor:
            replicate_counts[condition] += 1

        # DataSpec historically types paths as strings.  Preserve simple
        # single-pair paths here; lane lists remain in the raw artifact
        # manifest consumed by quantification.
        r1 = merged.get("R1", merged.get("r1"))
        r2 = merged.get("R2", merged.get("r2"))
        if isinstance(r1, str) and isinstance(r2, str):
            paths[str(name)] = {"R1": r1, "R2": r2}

    return {
        "paths": paths,
        "total_samples": len(raw_samples),
        "replicates": dict(replicate_counts),
        "batches": sorted(set(batches)),
        "conditions": list(dict.fromkeys(conditions)),
    }


class IntentParser:
    """Translate a natural-language query into a conservative StudyManifest."""

    @staticmethod
    def parse_prompt_to_manifest(
        user_prompt: str,
        study_id: str = "AD_mouse_001",
        raw_artifact_uri: Optional[str] = None,
        sample_manifest: Optional[Mapping[str, Any]] = None,
    ) -> StudyManifest:
        prompt = str(user_prompt or "")
        prompt_lower = prompt.lower()

        # Prefer the specific anatomical phrase before generic ``cortex`` or
        # ``brain`` matching.  Adrenal cortex is not brain cortex and must not
        # acquire an Alzheimer label just because it contains ``ad``.
        if _has_any_phrase(prompt_lower, ["adrenal cortex", "adrenal_cortex", "肾上腺皮质"]):
            tissue = "adrenal_cortex"
        elif _has_any_phrase(prompt_lower, ["brain", "脑"]):
            tissue = "brain"
        elif _has_any_phrase(prompt_lower, ["hippocampus", "海马"]):
            tissue = "hippocampus"
        elif _has_any_token(prompt_lower, ["cortex"]) or "皮质" in prompt_lower:
            tissue = "cortex"
        elif _has_any_phrase(prompt_lower, ["adrenal", "肾上腺"]):
            tissue = "adrenal"
        else:
            tissue = "unknown"

        if _has_any_phrase(prompt_lower, ["human", "人类", "人源"]):
            species = "homo_sapiens"
        elif _has_any_phrase(prompt_lower, ["mouse", "mice", "murine", "小鼠"]):
            species = "mus_musculus"
        else:
            species = "unknown"

        has_ko = _has_ko(prompt_lower)
        is_kat8 = _is_kat8(prompt_lower)
        ad_context = _is_ad(prompt_lower)
        if is_kat8:
            disease = "Kat8_cKO_Developmental_Disruption"
            conditions = ["cKO", "con"]
        elif ad_context:
            disease = "Alzheimer"
            conditions = ["AD", "control"]
        elif _has_any_phrase(prompt_lower, ["healthy", "normal", "健康"]):
            disease = "healthy"
            conditions = []
        else:
            disease = None
            conditions = []

        target_cells: List[str] = []
        if _has_any_phrase(prompt_lower, ["microglia", "小胶质"]):
            target_cells.append("Microglia")
        if _has_any_phrase(prompt_lower, ["astrocyte", "星形胶质"]):
            target_cells.append("Astrocytes")
        if _has_any_phrase(prompt_lower, ["neuron", "神经元"]):
            target_cells.append("Neurons")
        if _has_any_phrase(prompt_lower, ["progenitor", "前体", "干细胞", "stem"]):
            target_cells.append("Progenitors")

        modalities = ["scRNA"]
        has_fastq = bool(
            re.search(r"\b(?:fastq|fq)(?:\.gz)?\b", prompt_lower)
            or _has_any_phrase(
                prompt_lower,
                ["cleandata", "双端测序", "测序数据", "原始测序", "raw reads"],
            )
        )
        if has_fastq:
            modalities.append("FASTQ")

        has_spatial = _has_any_phrase(
            prompt_lower,
            ["spatial", "空间", "visium", "stereoseq", "merfish"],
        )
        if has_spatial:
            modalities.append("spatial")
        has_perturbation = has_ko or is_kat8 or _has_any_phrase(
            prompt_lower,
            ["perturb", "crispr", "knockout", "敲除", "cko"],
        )
        if has_perturbation:
            modalities.append("perturbation")
        if _has_any_phrase(
            prompt_lower,
            ["communication", "cci", "cellchat", "细胞通讯", "配体受体", "ligand-receptor"],
        ):
            modalities.append("communication")

        is_dam = _has_any_phrase(
            prompt_lower,
            ["dam假说", "prior", "trem2-apoe", "trem2", "apoe"],
        )
        is_prior_guided = is_kat8 or is_dam
        user_hypotheses: List[str] = []
        if is_kat8:
            user_hypotheses.append(
                "Kat8 (Mof) H4K16ac epigenetic loss impairs cell cycle progression and lineage differentiation"
            )
        elif is_dam:
            user_hypotheses.append("DAM subpopulation regulation via TREM2-APOE axis")

        sample_info = _extract_sample_manifest(sample_manifest)
        default_raw_uri = (
            f"fastq://{study_id}/raw_reads/v1"
            if has_fastq
            else f"adata://{study_id}/raw/v1"
        )
        biological_unit = (
            "mouse" if species == "mus_musculus" else "donor" if species == "homo_sapiens" else "unknown"
        )

        manifest = StudyManifest(
            study_id=study_id,
            title=f"Single-Cell Study: {disease or 'unknown context'} in {species} {tissue}",
            biological_design=BiologicalDesign(
                species=species,
                tissue=tissue,
                disease=disease,
                conditions=conditions or sample_info["conditions"],
                target_cell_types=target_cells,
            ),
            experimental_design=ExperimentalDesign(
                biological_unit=biological_unit,
                batches=sample_info["batches"],
                total_samples=sample_info["total_samples"],
                donor_replicates_per_condition=sample_info["replicates"],
            ),
            data=DataSpec(
                modalities=modalities,
                raw_artifact_uri=raw_artifact_uri or default_raw_uri,
                has_raw_fastq=has_fastq,
                fastq_paths=sample_info["paths"],
                has_spatial_coordinates=has_spatial,
                has_rna_velocity=False,
            ),
            hypotheses=Hypotheses(user_provided=user_hypotheses),
            analysis_policy=AnalysisPolicy(
                discovery_mode=not is_prior_guided,
                prior_guided_analysis=is_prior_guided,
                strict_reproducibility=True,
                prefer_pseudobulk=True,
            ),
            reproducibility=ReproducibilityConfig(random_seed=42),
        )
        return manifest
