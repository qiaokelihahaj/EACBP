"""Create a small SYNTHETIC ONLY run for manual WebUI acceptance checks.

Run from the repository root with its scientific Python environment. All data,
configuration and results stay under outputs/webui_optimization_validation.
The generated counts exercise real pipeline methods, but carry no biological
interpretation. Each invocation creates a separate run through the WebUI API.
"""
from __future__ import annotations

import json
from pathlib import Path
import time

import anndata as ad
import numpy as np
import pandas as pd

from eacbp.webui.service import WebService


def main() -> int:
    workspace = Path(__file__).resolve().parents[1]
    output = workspace / "outputs" / "webui_optimization_validation"
    output.mkdir(parents=True, exist_ok=True)
    notice = "SYNTHETIC ONLY - software acceptance data; no biological interpretation."
    (output / "README.txt").write_text(notice + "\n", encoding="utf-8")
    source = output / "SYNTHETIC_ONLY_webui_acceptance.h5ad"

    rng = np.random.default_rng(20260921)
    n_cells, n_genes = 192, 160
    donors = np.repeat([f"synthetic_donor_{i + 1}" for i in range(8)], 24)
    conditions = np.repeat(["control", "treated"], n_cells // 2)
    cell_types = np.tile(np.repeat(["Synthetic_A", "Synthetic_B"], 12), 8)
    means = np.full((n_cells, n_genes), 3.0)
    means[cell_types == "Synthetic_A", 4:20] = 12.0
    means[cell_types == "Synthetic_B", 20:36] = 12.0
    means[conditions == "treated", 40:48] = 10.0
    means[:, :4] = 0.15
    means *= np.repeat(rng.uniform(0.85, 1.15, 8), 24)[:, None]
    counts = rng.poisson(means).astype(np.float32)
    genes = [f"MT-SYNTHETIC_{i}" if i < 4 else f"SYNTHETIC_GENE_{i:03d}"
             for i in range(n_genes)]
    data = ad.AnnData(
        X=counts,
        obs=pd.DataFrame({"condition": conditions, "donor": donors,
                          "cell_type": cell_types, "validation_only": True},
                         index=[f"synthetic_cell_{i:03d}" for i in range(n_cells)]),
        var=pd.DataFrame({"gene_name": genes}, index=genes),
    )
    data.layers["counts"] = counts.copy()
    data.uns["validation_notice"] = notice
    data.uns["is_simulated"] = True
    data.uns["data_origin"] = "synthetic"
    data.write_h5ad(source)

    form = {
        "data": str(source), "study_id": "SYNTHETIC_ONLY_webui_acceptance",
        "title": notice, "species": "synthetic_test", "tissue": "synthetic_test",
        "method_profile": "baseline", "condition_col": "condition",
        "donor_col": "donor", "condition_a": "treated", "condition_b": "control",
        "min_genes": 40, "max_mito_pct": 20.0,
    }
    (output / "SYNTHETIC_ONLY_settings.json").write_text(
        json.dumps({"webui_schema": 1, "form": form}, indent=2), encoding="utf-8")
    service = WebService(workspace, output / "runs")
    started = service.start(form)
    deadline = time.monotonic() + 180
    while time.monotonic() < deadline:
        job = next(row for row in service.manager.jobs() if row["id"] == started["job"]["id"])
        if job["status"] not in {"queued", "running"}:
            evidence = {"notice": notice, "source": str(source),
                        "run_id": started["run_id"], "job": job}
            (output / "acceptance_run.json").write_text(
                json.dumps(evidence, indent=2, ensure_ascii=False), encoding="utf-8")
            print(json.dumps({"notice": notice, "source": str(source),
                              "run_id": started["run_id"], "status": job["status"],
                              "run_dir": job["run_dir"]}, indent=2))
            return 0 if job["status"] == "success" else 1
        time.sleep(0.2)
    raise TimeoutError("Synthetic acceptance worker did not finish within 180 seconds")


if __name__ == "__main__":
    raise SystemExit(main())
