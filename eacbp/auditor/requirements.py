"""Require evidence that each contracted audit actually ran.

Aliases describe the same check, never an unrelated successful check. Existing
check severities still decide whether a measured warning rejects a result.
Missing required checks always reject it, independently of those severities.
"""
from eacbp.auditor.base import ValidationCheck, ValidationSeverity
import numpy as np
import pandas as pd
from eacbp.capabilities.sc_data import SCData
from eacbp.schemas.artifact import ArtifactType

ALIASES = {
    "finite_expression_check": "expression_finite_values",
    "silhouette_check": "clustering_separation_silhouette",
    "marker_coherence_check": "canonical_marker_representation",
    "subset_non_empty": "matrix_non_empty",
    "fdr_correction_check": "multiple_testing_correction",
}


def missing_required_checks(requirements, checks):
    names = {check.check_name for check in checks}
    missing = []
    for requested in dict.fromkeys(requirements):
        actual = ALIASES.get(requested, requested)
        if actual not in names:
            missing.append(ValidationCheck(
                check_name=f"required_check_missing:{requested}", passed=False,
                severity=ValidationSeverity.ERROR,
                message=f"Required audit '{requested}' was not executed (canonical ID '{actual}').",
            ))
    return missing


def run_payload_requirements(contract, result, registry, checks):
    """Execute structural and design assessments previously only declared."""
    reports = []
    outputs = [registry.get(uri) for uri in result.output_artifacts]
    inputs = [registry.get(uri) for uri in contract.input_artifacts]
    def as_data(items):
        return [p if isinstance(p, SCData) else SCData.from_dict(p) for m, p in items
                if m.type in (ArtifactType.ANNDATA, ArtifactType.SPATIAL_DATA)]
    out_data, in_data = as_data(outputs), as_data(inputs)
    tables = [p for m, p in outputs if isinstance(p, pd.DataFrame)]
    present = {c.check_name for c in checks}
    for name in dict.fromkeys(contract.validation_requirements):
        if ALIASES.get(name, name) in present:
            continue
        ok, message, metrics = None, "", {}
        if name == "fdr_correction_check" or name == "multiple_testing_correction":
            ok = bool(tables) and all("fdr_q_value" in table for table in tables)
            if ok:
                for table in tables:
                    q = pd.to_numeric(table["fdr_q_value"], errors="coerce")
                    ok = ok and q.between(0, 1).all()
                    if "p_value" in table and not table.empty:
                        from scipy.stats import false_discovery_control
                        p = pd.to_numeric(table["p_value"], errors="coerce")
                        ok = ok and p.between(0, 1).all()
                        if p.between(0, 1).all():
                            ok = ok and np.allclose(q, false_discovery_control(p.to_numpy()), atol=1e-6)
            name = "multiple_testing_correction"
            message = "Checked adjusted significance range and recomputed BH where raw p-values are available."
        elif name == "silhouette_check" and contract.capability == "spatial_domain":
            ok = bool(out_data) and "X_spatial_pca" in out_data[0].obsm and "spatial_domain" in out_data[0].obs
            if ok:
                from sklearn.metrics import silhouette_score
                d = out_data[0]
                labels = d.obs["spatial_domain"].astype(str)
                score = float(silhouette_score(d.obsm["X_spatial_pca"], labels, sample_size=min(200, d.n_obs), random_state=42)) if 1 < labels.nunique() < d.n_obs else 0.0
                metrics = {"silhouette_score": score}
                ok = np.isfinite(score)
            name = "clustering_separation_silhouette"
            message = "Recomputed spatial-domain silhouette; separation is descriptive, not biological validation."
        elif name == "sample_count_check":
            ok = bool(in_data) and all(d.n_obs > 0 for d in in_data)
            if out_data and in_data:
                ok = ok and all(d.n_obs == in_data[0].n_obs for d in out_data)
            if contract.capability == "dataset_audit" and in_data and tables:
                ok = ok and len(tables[0]) == 1 and "n_cells" in tables[0] and int(tables[0].iloc[0]["n_cells"]) == in_data[0].n_obs
            message = "Compared recorded sample/cell dimensions against input artifacts."
        elif name == "replicate_check" and contract.capability == "dataset_audit":
            ok = bool(in_data) and bool(tables) and "min_replicates_per_condition" in tables[0]
            if ok:
                obs = in_data[0].obs
                condition = contract.parameters.get("condition_col", "condition")
                donor = contract.parameters.get("donor_col") or next((c for c in ("donor_id", "donor", "mouse_id", "sample_id", "sample") if c in obs), None)
                counts = obs.groupby(condition, observed=True)[donor].nunique() if condition in obs and donor else pd.Series(dtype=int)
                actual = int(counts.min()) if len(counts) else 0
                ok = int(tables[0].iloc[0]["min_replicates_per_condition"]) == actual
                metrics = {"min_replicates": actual, "replication_assessed": True}
            message = "Recomputed donor counts; an assessment does not establish sufficient replication."
        elif name == "retention_rate_check":
            ok = bool(out_data and in_data) and 0 < out_data[0].n_obs <= in_data[0].n_obs
            if ok:
                metrics = {"retention_rate": out_data[0].n_obs / in_data[0].n_obs}
            message = "Checked retained cells against the original input cell count."
        elif name == "hvg_count_check":
            ok = bool(out_data) and "highly_variable" in out_data[0].var
            if ok:
                flags = out_data[0].var["highly_variable"]
                ok = not flags.isna().any() and flags.isin([True, False]).all() and 0 < int(flags.sum()) <= out_data[0].n_vars
            message = "Checked the output highly-variable-gene mask."
        elif name == "finite_embedding_check":
            ok = bool(out_data) and all(bool(d.obsm) and all(np.isfinite(e.data if hasattr(e, "tocsr") else np.asarray(e)).all() and e.shape[0] == d.n_obs for e in d.obsm.values()) for d in out_data)
            message = "Checked all output embedding values and cell dimensions."
        elif name == "batch_mixing_audit":
            ok = bool(out_data) and "X_pca" in out_data[0].obsm
            if ok:
                from scipy.spatial import cKDTree
                d = out_data[0]
                embedding = np.asarray(d.obsm["X_pca"])
                ok = np.isfinite(embedding).all() and len(embedding) > 0
                if ok:
                    labels = d.obs["batch"].astype(str).to_numpy() if "batch" in d.obs else np.repeat("unknown", d.n_obs)
                    k = min(11, d.n_obs)
                    ix = cKDTree(embedding).query(embedding[:1000], k=k)[1]
                    score = float((labels[ix[:, 1:]] != labels[:len(ix), None]).mean()) if k > 1 else 0.0
                    metrics = {"neighbor_cross_batch_fraction": score, "batch_labels_available": "batch" in d.obs}
            message = "Recomputed neighbor batch mixing; this does not prove batch effects were removed."
        elif name == "dynamic_genes_check":
            ok = bool(tables) and {"gene", "fdr_q_value"}.issubset(tables[0].columns)
            if ok:
                q = pd.to_numeric(tables[0]["fdr_q_value"], errors="coerce")
                ok = q.between(0, 1).all()
            message = "Checked dynamic-gene identifiers and adjusted significance, allowing an empty result."
        elif name == "state_transition_stochasticity_check":
            matrices = [np.asarray(v) for d in out_data for k, v in d.obsp.items() if "transition" in k]
            matrices.extend(np.asarray(value["transition_matrix"]) for d in out_data for value in d.uns.values() if isinstance(value, dict) and "transition_matrix" in value)
            ok = bool(matrices) and all(np.isfinite(a).all() and (a >= 0).all() and np.allclose(a.sum(axis=1), 1, atol=1e-6) for a in matrices)
            message = "Checked stored transition matrices for finite nonnegative row-stochastic values."
        if ok is not None:
            reports.append(ValidationCheck(check_name=name, passed=bool(ok), severity=ValidationSeverity.ERROR,
                                           message=message, metrics=metrics))
    return reports
