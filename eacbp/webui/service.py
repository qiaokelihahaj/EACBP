"""Form validation and read-only views for a single-user local interface."""
from __future__ import annotations

from copy import deepcopy
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

from pydantic import BaseModel, ConfigDict, Field, StrictBool

from .jobs import JobManager, contained


class StudyForm(BaseModel):
    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)
    data: str = Field(min_length=1)
    study_id: str = Field(min_length=1, max_length=80, pattern=r"^[A-Za-z0-9][A-Za-z0-9_-]*$")
    title: str = ""
    species: str = Field(min_length=1)
    tissue: str = Field(min_length=1)
    target_cell_types: list[str] = Field(default_factory=list)
    method_profile: str = "standard"
    advanced_analysis: StrictBool = False
    condition_a: str = ""
    condition_b: str = ""
    condition_col: str = ""
    donor_col: str = ""
    batch_col: str = ""
    paired: StrictBool = False
    root_cell_id: str = ""
    min_genes: int = Field(default=100, ge=0)
    max_mito_pct: float = Field(default=20, ge=0, le=100)
    config_overrides: dict = Field(default_factory=dict)


def merge(base: dict, overrides: dict):
    result = deepcopy(base)
    for key, value in overrides.items():
        if isinstance(value, dict) and isinstance(result.get(key), dict):
            result[key] = merge(result[key], value)
        else:
            result[key] = deepcopy(value)
    return result


class WebService:
    def __init__(self, workspace: Path, runs_dir: Path):
        self.manager = JobManager(workspace, runs_dir)

    def settings(self):
        packages = {}
        for name in ("anndata", "scanpy", "harmonypy", "leidenalg", "cellrank", "pydeseq2",
                     "decoupler", "celltypist", "liana"):
            try:
                packages[name] = version(name)
            except PackageNotFoundError:
                packages[name] = None
        return {"workspace": str(self.manager.workspace), "runs_dir": str(self.manager.runs_dir),
                "packages": packages, "max_concurrent_jobs": 1}

    def browse(self, relative="."):
        root = self.manager.workspace
        path = contained(root, relative)
        if not path.is_dir():
            raise ValueError("请选择文件夹")
        entries = []
        # Hide generated/scratch trees by default, while allowing an explicit path.
        for child in path.iterdir():
            if child.name.startswith((".", "temp_", "tmp_", "__pycache__")):
                continue
            try:
                resolved = contained(root, child)
                if child.is_dir() or child.suffix.casefold() == ".h5ad":
                    entries.append({"name": child.name, "path": str(resolved),
                                    "directory": child.is_dir(),
                                    "size": child.stat().st_size if child.is_file() else None})
            except (OSError, ValueError):
                continue
        entries.sort(key=lambda row: (not row["directory"], row["name"].casefold()))
        return {"path": str(path), "parent": str(path.parent) if path != root else None,
                "entries": entries[:300], "truncated": len(entries) > 300}

    def data_path(self, value):
        path = contained(self.manager.workspace, value)
        if path.suffix.casefold() != ".h5ad" or not path.is_file():
            raise ValueError("请选择现有的 .h5ad 文件（路径相对于工作目录）")
        return path

    def dataset(self, value):
        import anndata as ad
        path = self.data_path(value)
        data = ad.read_h5ad(path, backed="r")
        try:
            columns = {}
            for name in data.obs.columns:
                series = data.obs[name]
                unique = series.dropna().astype(str).unique()
                columns[str(name)] = {"unique": int(len(unique)), "values": unique[:12].tolist(),
                                      "missing": int(series.isna().sum())}
            return {"path": str(path), "n_cells": data.n_obs, "n_genes": data.n_vars,
                    "columns": columns, "layers": list(data.layers.keys()),
                    "spatial": "spatial" in data.obsm,
                    "notice": "只读取结构与元数据；尚未验证计数、实验设计或分析适用性。"}
        finally:
            data.file.close()

    def prepare(self, payload):
        from eacbp.cli import _normalise_config_paths
        from eacbp.schemas.runtime import RunConfig
        from eacbp.schemas.study import StudyManifest

        form = StudyForm.model_validate(payload)
        data = self.data_path(form.data)
        if bool(form.condition_a) != bool(form.condition_b):
            raise ValueError("对比条件 A 和 B 必须同时填写，或同时留空")
        if form.condition_a and form.condition_a == form.condition_b:
            raise ValueError("对比条件 A 和 B 不能相同")
        if any(not label.strip() for label in form.target_cell_types):
            raise ValueError("目标细胞类型不能包含空值")
        params = {"qc": {"min_genes": form.min_genes, "max_mito_pct": form.max_mito_pct}}
        for capability in ("dataset_audit", "deg", "differential_abundance"):
            params[capability] = {key: getattr(form, key) for key in
                                  ("condition_col", "donor_col", "batch_col") if getattr(form, key)}
        if form.batch_col:
            params["integration"] = {"batch_col": form.batch_col}
        if form.condition_a:
            for capability in ("deg", "differential_abundance"):
                params[capability].update(condition_a=form.condition_a, condition_b=form.condition_b)
        if form.paired:
            params["deg"]["paired"] = True
        if form.root_cell_id:
            params["trajectory_inference"] = {"root_cell_id": form.root_cell_id}
        config = merge({"mode": "real", "method_profile": form.method_profile,
                        "advanced_analysis": form.advanced_analysis, "capability_parameters": params},
                       form.config_overrides)
        config = RunConfig.from_mapping(_normalise_config_paths(config, base_dir=self.manager.workspace)).to_mapping()
        if config.get("mode") != "real" or config.get("resume") or config.get("has_fastq"):
            raise ValueError("新建界面仅接收真实 h5ad 输入；恢复请使用历史运行按钮")
        selected_deg = config.get("method_overrides", {}).get("deg") or (
            "pydeseq2_pseudobulk_v1" if config.get("advanced_analysis") else "donor_pseudobulk_welch_v1"
        )
        if config.get("capability_parameters", {}).get("deg", {}).get("paired") and selected_deg not in {
            "pydeseq2_pseudobulk_v1", "pydeseq2_deg_v1", "pydeseq2_donor_pseudobulk_v1",
        }:
            raise ValueError("配对设计需要开启 PyDESeq2 高级分析")
        manifest = StudyManifest.model_validate({
            "study_id": form.study_id, "title": form.title or form.study_id,
            "biological_design": {"species": form.species, "tissue": form.tissue,
                                  "target_cell_types": form.target_cell_types},
            "experimental_design": {"biological_unit": "donor", "total_samples": 0},
            "data": {"raw_artifact_uri": f"adata://{form.study_id}/raw/v1"},
        })
        return form, data, manifest, config

    def preview(self, payload):
        from eacbp.orchestrator.planning import preview_study_plan
        form, data, manifest, config = self.prepare(payload)
        plan = preview_study_plan(manifest, config)
        required = {"anndata"}
        by_method = {"harmonypy_v1": ["harmonypy", "scanpy"],
                     "scanpy_leiden_umap_v1": ["scanpy", "leidenalg"],
                     "scanpy_dpt_v1": ["scanpy"], "cellrank_fate_v1": ["cellrank"],
                     "pydeseq2_pseudobulk_v1": ["pydeseq2"],
                     "pydeseq2_leave_one_donor_out_v1": ["pydeseq2"],
                     "decoupler_ulm_v2": ["decoupler"]}
        for task in plan["tasks"]:
            required.update(by_method.get(task.get("method"), []))
        installed = self.settings()["packages"]
        missing = sorted(package for package in required if not installed.get(package))
        return {"form": form.model_dump(), "data": str(data),
                "manifest": manifest.model_dump(mode="json"), "config": config,
                "plan": plan, "missing_packages": missing,
                "can_run": plan["valid"] and not missing}

    def start(self, payload):
        prepared = self.preview(payload)
        if not prepared["can_run"]:
            raise ValueError("计划包含错误或缺少依赖，请先修复预览中的问题")
        job = self.manager.submit(
            "run", run_dir=self.manager.new_run_dir(prepared["manifest"]["study_id"]),
            manifest=prepared["manifest"], config=prepared["config"], data=prepared["data"],
        )
        return {"job": job, "run_id": Path(job["run_dir"]).relative_to(self.manager.runs_dir).as_posix()}
