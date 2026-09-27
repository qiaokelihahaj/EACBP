"""Form validation and read-only views for a single-user local interface."""
from __future__ import annotations

from copy import deepcopy
from contextlib import ExitStack
import csv
import hashlib
import json
from importlib.metadata import PackageNotFoundError, version
import os
from pathlib import Path
import secrets
import sys
import tempfile
import threading
import time
from typing import Literal
from uuid import uuid4

from pydantic import BaseModel, ConfigDict, Field, StrictBool

from .jobs import JobManager, contained, lease, read_json, write_json
from .preferences import DirectorySettings, absolute_directory


MAX_IMPORT_FILE_BYTES = 64_000_000
MAX_IMPORT_MATRIX_VALUES = 2_000_000
MAX_IMPORT_CELLS = 250_000
MAX_IMPORT_GENES = 100_000


class StudyForm(BaseModel):
    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)
    data: str = ""
    study_id: str = Field(default="study_001", min_length=1, max_length=80, pattern=r"^[A-Za-z0-9][A-Za-z0-9_-]*$")
    title: str = ""
    species: str = ""
    tissue: str = ""
    research_purpose: Literal["overview", "deg", "abundance", "trajectory", "batch", "annotation"] = "overview"
    target_cell_types: list[str] = Field(default_factory=list)
    method_profile: str = "standard"
    advanced_analysis: StrictBool = False
    condition_a: str = ""
    condition_b: str = ""
    condition_col: str = ""
    sample_col: str = ""
    donor_col: str = ""
    batch_col: str = ""
    paired: StrictBool = False
    root_cell_id: str = ""
    design_formula: str = ""
    counts_layer: str = "counts"
    allow_x_as_counts: StrictBool = False
    min_donors: int = Field(default=2, ge=2)
    alpha: float = Field(default=0.05, gt=0, lt=1)
    celltypist_model_path: str = ""
    annotation_use_as_cell_type: StrictBool = False
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
    def __init__(self, workspace: Path, runs_dir: Path, *, settings_file: Path | None = None):
        self.directory_guard = threading.RLock()
        self.settings_file = Path(settings_file).resolve() if settings_file is not None else None
        self.manager = JobManager(workspace, runs_dir)
        self._import_previews = {}

    def settings(self):
        packages = {}
        for name in ("anndata", "scanpy", "harmonypy", "leidenalg", "cellrank", "pydeseq2",
                     "decoupler", "celltypist", "liana"):
            try:
                packages[name] = version(name)
            except PackageNotFoundError:
                packages[name] = None
        return {"workspace": str(self.manager.workspace), "runs_dir": str(self.manager.runs_dir),
                "package_dir": str(Path(__file__).resolve().parents[1]), "python_executable": sys.executable,
                "settings_file": str(self.settings_file) if self.settings_file else None,
                "packages": packages, "max_concurrent_jobs": 1}

    def update_directories(self, payload):
        """Publish both roots together, only while old and new job queues are idle."""
        form = DirectorySettings.model_validate(payload)
        workspace, runs_dir = absolute_directory(form.workspace), absolute_directory(form.runs_dir)
        if not workspace.is_dir():
            raise ValueError("输入数据目录不存在，请先创建或选择现有目录")
        # Check directory listing access without loading any input files.
        with os.scandir(workspace) as entries:
            next(entries, None)
        with self.directory_guard, ExitStack() as locks:
            current = self.manager
            locks.enter_context(current._guard)
            locks.enter_context(lease(current.jobs_dir.parent / "submit.lock"))
            if any(job["status"] in {"queued", "running", "unknown"} for job in current.jobs()):
                raise BlockingIOError("当前目录有运行中或状态待检查的任务，暂不能切换目录")
            candidate = JobManager(workspace, runs_dir)
            if runs_dir != current.runs_dir:
                locks.enter_context(lease(candidate.jobs_dir.parent / "submit.lock"))
                if any(job["status"] in {"queued", "running", "unknown"} for job in candidate.jobs()):
                    raise BlockingIOError("目标结果目录有运行中或状态待检查的任务，暂不能切换")
            with tempfile.TemporaryFile(dir=runs_dir) as probe:
                probe.write(b"eacbp directory access check")
                probe.flush()
            if self.settings_file is not None:
                write_json(self.settings_file, {"schema_version": 1, "workspace": str(workspace), "runs_dir": str(runs_dir)})
            # A failed save must leave the current roots untouched.
            self.manager = candidate
            self._import_previews.clear()
        return {"workspace": str(workspace), "runs_dir": str(runs_dir),
                "saved": self.settings_file is not None}

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
                if child.is_dir() or child.suffix.casefold() in {".h5ad", ".csv", ".tsv"}:
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
                columns[str(name)] = {"unique": int(len(unique)), "values": unique[:100].tolist(),
                                      "truncated": len(unique) > 100,
                                      "missing": int(series.isna().sum())}
            return {"path": str(path), "n_cells": data.n_obs, "n_genes": data.n_vars,
                    "columns": columns, "layers": list(data.layers.keys()),
                    "spatial": "spatial" in data.obsm,
                    "notice": "只读取结构与元数据；尚未验证计数、实验设计或分析适用性。"}
        finally:
            data.file.close()

    @staticmethod
    def _csv_digest(path: Path):
        digest = hashlib.sha256()
        with path.open("rb") as stream:
            while chunk := stream.read(1024 * 1024):
                digest.update(chunk)
        return digest.hexdigest()

    def _import_file(self, value, role):
        path = contained(self.manager.workspace, value)
        if path.suffix.casefold() not in {".csv", ".tsv"} or not path.is_file():
            raise ValueError(f"{role} 必须是工作目录内现有的 CSV 或 TSV 文件。")
        if path.stat().st_size > MAX_IMPORT_FILE_BYTES:
            raise ValueError(f"{role} 超过 64 MB 导入预算；请先整理或拆分输入文件。")
        return path

    @staticmethod
    def _csv_shape(path: Path, *, role: str, max_rows: int, max_columns: int,
                   max_values: int | None = None, count_missing=False):
        delimiter = "\t" if path.suffix.casefold() == ".tsv" else ","
        try:
            with path.open("r", encoding="utf-8-sig", newline="") as stream:
                reader = csv.reader(stream, delimiter=delimiter, strict=True)
                try:
                    header = next(reader)
                except StopIteration as exc:
                    raise ValueError(f"{role} 文件为空。") from exc
                if len(header) < 2 or len(header) > max_columns:
                    raise ValueError(f"{role} 表头需包含 ID 列和数据列，且列数不超过 {max_columns}。")
                if any(not name or not name.strip() for name in header):
                    raise ValueError(f"{role} 表头包含空列名；请明确填写唯一的列名。")
                seen_headers, duplicate_headers = set(), set()
                for name in header:
                    if name in seen_headers:
                        duplicate_headers.add(name)
                    seen_headers.add(name)
                duplicates = sorted(duplicate_headers)
                if duplicates:
                    raise ValueError(f"{role} 表头有重复列名：{duplicates[:5]}。")
                identifiers, rows = [], 0
                missing = {column: 0 for column in header[1:]} if count_missing else {}
                seen = set()
                for row in reader:
                    rows += 1
                    if rows > max_rows:
                        raise ValueError(f"{role} 行数超过导入预算（最多 {max_rows}）。")
                    if len(row) != len(header):
                        raise ValueError(f"{role} 第 {reader.line_num} 行有 {len(row)} 列，应为 {len(header)} 列。")
                    identifier = row[0]
                    if not identifier or not identifier.strip():
                        raise ValueError(f"{role} 第 {reader.line_num} 行 ID 缺失。")
                    if identifier in seen:
                        raise ValueError(f"{role} 包含重复 ID：{identifier[:120]}。")
                    seen.add(identifier)
                    identifiers.append(identifier)
                    if count_missing:
                        for column, value in zip(header[1:], row[1:]):
                            if value == "":
                                missing[column] += 1
                    if max_values is not None and rows * (len(header) - 1) > max_values:
                        raise ValueError(f"计数矩阵超过 {MAX_IMPORT_MATRIX_VALUES:,} 个数据值的内存预算。")
                if rows == 0:
                    raise ValueError(f"{role} 文件没有数据行。")
                return {"delimiter": delimiter, "header": header, "ids": identifiers,
                        "rows": rows, "missing": missing}
        except UnicodeDecodeError as exc:
            raise ValueError(f"{role} 必须使用 UTF-8 或 UTF-8 BOM 编码。") from exc
        except csv.Error as exc:
            raise ValueError(f"{role} CSV/TSV 格式无效：{exc}。") from exc

    def _validate_import_sources(self, payload):
        if not isinstance(payload, dict):
            raise ValueError("导入预检需要计数矩阵、元数据路径和明确矩阵方向。")
        orientation = payload.get("orientation")
        if orientation not in {"cells_by_genes", "genes_by_cells"}:
            raise ValueError("请选择矩阵方向：细胞为行、基因为列，或基因为行、细胞为列。")
        counts = self._import_file(payload.get("counts", ""), "计数矩阵")
        metadata = self._import_file(payload.get("metadata", ""), "细胞元数据")
        matrix_cells_as_rows = orientation == "cells_by_genes"
        matrix = self._csv_shape(counts, role="计数矩阵",
                                 max_rows=MAX_IMPORT_CELLS if matrix_cells_as_rows else MAX_IMPORT_GENES,
                                 max_columns=(MAX_IMPORT_GENES if matrix_cells_as_rows else MAX_IMPORT_CELLS) + 1,
                                 max_values=MAX_IMPORT_MATRIX_VALUES)
        observations = self._csv_shape(metadata, role="细胞元数据", max_rows=MAX_IMPORT_CELLS,
                                       max_columns=501, count_missing=True)
        if orientation == "cells_by_genes":
            cell_ids, gene_ids = matrix["ids"], matrix["header"][1:]
        else:
            gene_ids, cell_ids = matrix["ids"], matrix["header"][1:]
        if len(cell_ids) > MAX_IMPORT_CELLS or len(gene_ids) > MAX_IMPORT_GENES:
            raise ValueError(f"矩阵最多支持 {MAX_IMPORT_CELLS:,} 个细胞和 {MAX_IMPORT_GENES:,} 个基因。")
        if len(observations["header"]) < 2:
            raise ValueError("细胞元数据至少需要一个 ID 列和一个元数据列。")
        count_set, metadata_set = set(cell_ids), set(observations["ids"])
        if count_set != metadata_set:
            missing = sorted(count_set - metadata_set)[:5]
            extra = sorted(metadata_set - count_set)[:5]
            raise ValueError("计数矩阵与元数据的细胞 ID 集合必须完全相同；"
                             f"元数据缺少 {len(count_set - metadata_set)} 个，"
                             f"另有 {len(metadata_set - count_set)} 个未匹配。"
                             f"缺少示例：{missing}；多余示例：{extra}。")
        if len(cell_ids) * len(gene_ids) > MAX_IMPORT_MATRIX_VALUES:
            raise ValueError(f"计数矩阵超过 {MAX_IMPORT_MATRIX_VALUES:,} 个数据值的内存预算。")
        return {"counts_path": counts, "metadata_path": metadata, "orientation": orientation,
                "counts_shape": matrix, "metadata_shape": observations,
                "cell_ids": cell_ids, "gene_ids": gene_ids,
                "counts_sha256": self._csv_digest(counts),
                "metadata_sha256": self._csv_digest(metadata)}

    @staticmethod
    def _load_import_data(inspected):
        import numpy as np
        import pandas as pd

        def read_csv(path):
            sep = "\t" if path.suffix.casefold() == ".tsv" else ","
            return pd.read_csv(path, sep=sep, encoding="utf-8-sig", dtype=str,
                               keep_default_na=False)

        counts_frame = read_csv(inspected["counts_path"])
        metadata_frame = read_csv(inspected["metadata_path"])
        if inspected["orientation"] == "cells_by_genes":
            cell_ids = counts_frame.iloc[:, 0].tolist()
            gene_ids = counts_frame.columns[1:].tolist()
            numeric = counts_frame.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")
        else:
            cell_ids = counts_frame.columns[1:].tolist()
            gene_ids = counts_frame.iloc[:, 0].tolist()
            numeric = counts_frame.iloc[:, 1:].apply(pd.to_numeric, errors="coerce").T
        if cell_ids != inspected["cell_ids"] or gene_ids != inspected["gene_ids"]:
            raise ValueError("計數矩陣在读取期间发生变化；请重新预检。")
        values = numeric.to_numpy(dtype=np.float64, copy=False)
        if (not np.isfinite(values).all() or (values < 0).any()
                or (values != np.floor(values)).any()
                or (values > np.iinfo(np.int32).max).any()):
            raise ValueError("计数矩阵含缺失、非有限、负数、非整数或超出 int32 的数值；请提供原始非负整数计数。")
        matrix = values.astype(np.int32, copy=False)
        id_column = metadata_frame.columns[0]
        if metadata_frame[id_column].tolist() != inspected["metadata_shape"]["ids"]:
            raise ValueError("细胞元数据在读取期间发生变化；请重新预检。")
        observations = metadata_frame.set_index(id_column, verify_integrity=True)
        observations = observations.loc[cell_ids].copy()
        observations = observations.mask(observations.eq(""), pd.NA)
        observations.index = pd.Index(cell_ids, dtype=object, name=id_column)
        return matrix, observations, cell_ids, gene_ids

    def preview_import(self, payload):
        import numpy as np

        inspected = self._validate_import_sources(payload)
        matrix, _, cell_ids, gene_ids = self._load_import_data(inspected)
        if (self._csv_digest(inspected["counts_path"]) != inspected["counts_sha256"]
                or self._csv_digest(inspected["metadata_path"]) != inspected["metadata_sha256"]):
            raise ValueError("预检过程中源文件发生变化；请稍后重试。")
        now = time.time()
        self._import_previews = {key: value for key, value in self._import_previews.items()
                                 if value["expires"] > now}
        if len(self._import_previews) >= 8:
            self._import_previews.pop(next(iter(self._import_previews)))
        token = secrets.token_urlsafe(24)
        self._import_previews[token] = {**inspected, "expires": now + 1800}
        n_cells, n_genes = len(inspected["cell_ids"]), len(inspected["gene_ids"])
        return {"preview_token": token, "counts_path": str(inspected["counts_path"]),
                "metadata_path": str(inspected["metadata_path"]),
                "orientation": inspected["orientation"], "n_cells": n_cells, "n_genes": n_genes,
                "metadata_columns": inspected["metadata_shape"]["header"][1:],
                "metadata_missing": {column: int(inspected["metadata_shape"]["missing"].get(column, 0))
                                     for column in inspected["metadata_shape"]["header"][1:]},
                "metadata_rows_reordered": inspected["cell_ids"] != inspected["metadata_shape"]["ids"],
                "counts_sha256": inspected["counts_sha256"],
                "metadata_sha256": inspected["metadata_sha256"],
                "count_check": "全部计数为有限、非负 int32 整数；元数据空单元格将转换为缺失值；输入未写入。",
                "count_range": {"minimum": int(matrix.min()), "maximum": int(matrix.max()),
                                "nonzero": int(np.count_nonzero(matrix))},
                "alignment": "数据和元数据细胞 ID 集合完全匹配；元数据按 ID 对齐。",
                "limits": {"each_file_bytes": MAX_IMPORT_FILE_BYTES,
                           "matrix_values": MAX_IMPORT_MATRIX_VALUES}}

    def confirm_import(self, payload):
        token = payload.get("preview_token") if isinstance(payload, dict) else None
        preview = self._import_previews.pop(token, None)
        if not preview or preview["expires"] <= time.time():
            raise ValueError("导入预览已过期或不存在；请重新预检两个源文件。")
        inspected = self._validate_import_sources({
            "counts": str(preview["counts_path"]), "metadata": str(preview["metadata_path"]),
            "orientation": preview["orientation"]})
        if (inspected["counts_sha256"] != preview["counts_sha256"]
                or inspected["metadata_sha256"] != preview["metadata_sha256"]):
            raise ValueError("预览后源文件发生变化；请重新预检并确认。")

        import anndata as ad
        import pandas as pd
        matrix, observations, cell_ids, gene_ids = self._load_import_data(inspected)
        if (self._csv_digest(inspected["counts_path"]) != preview["counts_sha256"]
                or self._csv_digest(inspected["metadata_path"]) != preview["metadata_sha256"]):
            raise ValueError("确认导入期间源文件发生变化；没有创建 h5ad，请重新预检。")
        variables = pd.DataFrame(index=pd.Index(gene_ids, dtype=object,
                                                  name=inspected["counts_shape"]["header"][0]
                                                  if preview["orientation"] == "genes_by_cells"
                                                  else None))
        output_dir = contained(self.manager.workspace, self.manager.workspace / ".eacbp" / "imports")
        output_dir.mkdir(parents=True, exist_ok=True)
        output_dir = contained(self.manager.workspace, output_dir)
        output = output_dir / f"{time.strftime('%Y%m%d-%H%M%S')}-{uuid4().hex[:10]}.h5ad"
        temporary = output.with_name(output.name + ".tmp")
        data = ad.AnnData(X=matrix, obs=observations, var=variables)
        data.layers["counts"] = matrix.copy()
        data.uns["eacbp_csv_import"] = {
            "schema_version": 1,
            "orientation": preview["orientation"],
            "counts_file_name": inspected["counts_path"].name,
            "metadata_file_name": inspected["metadata_path"].name,
            "counts_sha256": inspected["counts_sha256"],
            "metadata_sha256": inspected["metadata_sha256"],
            "cell_ids_matched_by": "exact ID set; metadata reordered to count matrix cell IDs",
            "count_validation": "finite nonnegative int32 integers",
        }
        try:
            data.write_h5ad(temporary)
            if output.exists():
                raise FileExistsError("生成的导入路径已存在；请再次确认导入以创建新文件。")
            os.replace(temporary, output)
        finally:
            temporary.unlink(missing_ok=True)
        return {"data": str(output), "n_cells": len(cell_ids), "n_genes": len(gene_ids),
                "notice": "已创建新 h5ad，未覆盖源文件。请检查结构并重新预览分析计划。"}

    def _configuration(self, form, *, validate_contrast=True):
        from eacbp.application.study_service import normalise_config_paths
        from eacbp.schemas.runtime import RunConfig

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
        if form.advanced_analysis:
            params["deg"].update(
                counts_layer=form.counts_layer,
                allow_x_as_counts=form.allow_x_as_counts,
                min_donors=form.min_donors,
                alpha=form.alpha,
            )
            if form.design_formula:
                params["deg"]["design_formula"] = form.design_formula
        base = {"mode": "real", "method_profile": form.method_profile,
                "advanced_analysis": form.advanced_analysis, "capability_parameters": params}
        if form.celltypist_model_path or form.annotation_use_as_cell_type:
            base["analysis_extensions"] = {"cell_annotation": {
                **({"model_path": form.celltypist_model_path} if form.celltypist_model_path else {}),
                "use_as_cell_type": form.annotation_use_as_cell_type,
            }}
        config = merge(base, form.config_overrides)
        config = RunConfig.from_mapping(normalise_config_paths(config, base_dir=self.manager.workspace)).to_mapping()
        if validate_contrast:
            for capability in ("deg", "differential_abundance"):
                values = config.get("capability_parameters", {}).get(capability, {}) or {}
                condition_a, condition_b = values.get("condition_a"), values.get("condition_b")
                if bool(condition_a) != bool(condition_b):
                    raise ValueError("对比条件 A 和 B 必须同时填写，或同时留空")
                if condition_a and condition_a == condition_b:
                    raise ValueError("对比条件 A 和 B 不能相同")
        return config

    def reuse_configuration(self, run_id):
        """Project a successful run's final config into explicit, editable form fields."""
        run = self.manager.run_path(run_id)
        if not run.is_dir():
            raise FileNotFoundError("运行记录不存在")
        latest = next((job for job in self.manager.jobs() if job.get("run_dir") == str(run)), None)
        if latest and latest["status"] in {"queued", "running", "unknown"}:
            raise BlockingIOError("该运行仍在执行或状态待检查，完成后再复用配置。")
        saved = read_json(run / "run_config.json", {})
        if saved.get("status") != "success" or saved.get("import_completed") is not True:
            raise ValueError("只有已成功完成数据导入的运行才能复用配置。")
        original_config = saved.get("config")
        manifest = saved.get("manifest")
        if not isinstance(original_config, dict) or not isinstance(manifest, dict):
            raise ValueError("历史运行缺少完整的配置或研究元数据，无法安全复用。")
        if original_config.get("mode", "real") != "real" or original_config.get("has_fastq"):
            raise ValueError("此历史运行不是可复用的真实 h5ad 配置。")

        original_config = deepcopy(original_config)
        preserved = deepcopy(original_config)
        for key in ("resume", "has_fastq", "data"):
            preserved.pop(key, None)

        def take(capabilities, key, default=None):
            for capability in capabilities:
                values = original_config.get("capability_parameters", {}).get(capability, {}) or {}
                if key in values and values[key] is not None:
                    return values[key]
            return default

        form = {
            "data": "",
            "study_id": "",
            "title": "",
            "species": "",
            "tissue": "",
            "research_purpose": "overview",
            "target_cell_types": [],
            "method_profile": original_config.get("method_profile", "standard"),
            "advanced_analysis": original_config.get("advanced_analysis", False),
            "condition_a": take(("deg", "differential_abundance"), "condition_a", ""),
            "condition_b": take(("deg", "differential_abundance"), "condition_b", ""),
            "condition_col": take(("deg", "differential_abundance", "dataset_audit"), "condition_col", ""),
            "sample_col": "",
            "donor_col": take(("deg", "differential_abundance", "dataset_audit"), "donor_col", ""),
            "batch_col": take(("deg", "differential_abundance", "dataset_audit", "integration"), "batch_col", ""),
            "paired": take(("deg",), "paired", False),
            "root_cell_id": take(("trajectory_inference",), "root_cell_id", ""),
            "counts_layer": take(("deg",), "counts_layer", "counts"),
            "allow_x_as_counts": take(("deg",), "allow_x_as_counts", False),
            "min_donors": take(("deg",), "min_donors", 2),
            "alpha": take(("deg",), "alpha", 0.05),
            "design_formula": take(("deg",), "design_formula", ""),
            "celltypist_model_path": "",
            "annotation_use_as_cell_type": False,
            "min_genes": take(("qc",), "min_genes", 100),
            "max_mito_pct": take(("qc",), "max_mito_pct", 20),
        }
        biological = manifest.get("biological_design", {}) or {}
        form["species"] = biological.get("species") or ""
        form["tissue"] = biological.get("tissue") or ""
        form["target_cell_types"] = list(biological.get("target_cell_types") or [])
        original_id = str(saved.get("study_id") or manifest.get("study_id") or "study_001")
        form["study_id"] = original_id[:73] + "_reuse"
        form["title"] = str(manifest.get("title") or original_id) + " (复用配置)"

        original_extensions = original_config.get("analysis_extensions") or {}
        if not isinstance(original_extensions, dict):
            original_extensions = {}
        cell_annotation = original_extensions.get("cell_annotation") or {}
        if not isinstance(cell_annotation, dict):
            cell_annotation = {}
        cleared_external_inputs = []
        model_path = cell_annotation.get("model_path")
        if model_path:
            if Path(str(model_path)).is_file():
                form["celltypist_model_path"] = str(model_path)
            else:
                cleared_external_inputs.append({"name": "CellTypist 模型路径", "path": str(model_path),
                                               "reason": "原配置引用的文件当前不存在，已从复用配置中移除。"})
                current_extensions = preserved.get("analysis_extensions") or {}
                if isinstance(current_extensions, dict):
                    current_annotation = current_extensions.get("cell_annotation") or {}
                    if isinstance(current_annotation, dict):
                        current_annotation.pop("model_path", None)
        research_context = read_json(run / "webui_context.json", {}) or {}
        purpose_saved = research_context.get("research_purpose") in {
            "overview", "deg", "abundance", "trajectory", "batch", "annotation"}
        if purpose_saved:
            form["research_purpose"] = research_context["research_purpose"]
        annotation_saved = cell_annotation
        if "use_as_cell_type" in annotation_saved:
            form["annotation_use_as_cell_type"] = bool(annotation_saved["use_as_cell_type"])
        form["config_overrides"] = preserved
        unsupported_manifest_sections = [
            "实验设计字段（生物学单位、批次数和重复数）会按新输入重新生成",
            "数据规格（包括旧 raw_artifact_uri、模态和 FASTQ 路径）不会迁移",
            "假设、运行约束、分析策略和复现性设置会使用 WebUI 新建运行的默认值",
        ]
        if (manifest.get("biological_design") or {}).get("disease"):
            unsupported_manifest_sections.append("biological_design.disease 未映射到新建表单")
        return {"form": form,
                "original": {"run_id": run_id, "study_id": original_id,
                             "source_path": saved.get("source_path") or (saved.get("source") or {}).get("path"),
                             "config": original_config, "manifest": manifest},
                "research_purpose_saved": purpose_saved,
                "cleared_external_inputs": cleared_external_inputs,
                "unsupported_manifest_sections": unsupported_manifest_sections,
                "notice": "复用了可映射的历史配置。原数据路径已清空，需选择当前输入并重新检查、预览和确认。"}

    def prepare(self, payload):
        from eacbp.schemas.study import StudyManifest

        form = StudyForm.model_validate(payload)
        data = self.data_path(form.data)
        if not form.species:
            raise ValueError("请填写物种")
        if not form.tissue:
            raise ValueError("请填写组织")
        config = self._configuration(form)
        if config.get("mode") != "real" or config.get("resume") or config.get("has_fastq"):
            raise ValueError("新建分析需要真实 h5ad 输入；CSV/TSV 请先通过导入向导转换，恢复请使用历史运行按钮")
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

    @staticmethod
    def _effective_column(parameters, key, fallback, *, empty_falls_back=False):
        if key not in parameters:
            return fallback
        value = parameters[key]
        if value is None or (empty_falls_back and value == ''):
            return fallback
        return str(value).strip()

    def design(self, payload, *, plan=None, prepared_config=None):
        """Summarize obs metadata and check selected design requirements."""
        import anndata as ad

        form = StudyForm.model_validate(payload)
        path = self.data_path(form.data)
        config = prepared_config if prepared_config is not None else self._configuration(form, validate_contrast=bool(plan))
        raw_parameters = config.get("capability_parameters") or {}
        tasks = {}
        tasks_by_capability = {}
        if plan:
            for task in plan.get("tasks", []):
                tasks.setdefault(task.get("capability"), task)
                tasks_by_capability.setdefault(task.get("capability"), []).append(task)

        def parameters(capability):
            task = tasks.get(capability)
            if task is not None:
                return dict(task.get("parameters") or {})
            return dict(raw_parameters.get(capability) or {})

        def all_parameters(capability):
            planned = tasks_by_capability.get(capability, [])
            return [dict(task.get("parameters") or {}) for task in planned] or [
                dict(raw_parameters.get(capability) or {})
            ]

        audit = parameters("dataset_audit")
        deg = parameters("deg")
        abundance = parameters("differential_abundance")
        integration = parameters("integration")
        trajectory = parameters("trajectory_inference")
        audit_raw = dict(raw_parameters.get("dataset_audit") or {})
        deg_raw = dict(raw_parameters.get("deg") or {})
        abundance_raw = dict(raw_parameters.get("differential_abundance") or {})
        integration_raw = dict(raw_parameters.get("integration") or {})

        data = ad.read_h5ad(path, backed="r")
        try:
            obs = data.obs
            names = {str(name) for name in obs.columns}
            audit_donor_fallback = next((name for name in ("donor_id", "donor", "mouse_id", "sample_id", "sample") if name in names), "")
            stats_donor_fallback = next((name for name in ("donor", "donor_id", "mouse_id", "sample_id", "sample") if name in names), "")
            effective = {
                "audit_condition": self._effective_column(audit, "condition_col", "condition"),
                "deg_condition": self._effective_column(deg, "condition_col", "condition"),
                "abundance_condition": self._effective_column(abundance, "condition_col", "condition"),
                "audit_donor": self._effective_column(audit, "donor_col", audit_donor_fallback, empty_falls_back=True),
                "deg_donor": self._effective_column(deg, "donor_col", stats_donor_fallback),
                "abundance_donor": self._effective_column(abundance, "donor_col", stats_donor_fallback),
                "audit_batch": self._effective_column(audit, "batch_col", "batch" if "batch" in names else "", empty_falls_back=True),
                "integration_batch": self._effective_column(integration, "batch_col", "batch" if "batch" in names else "", empty_falls_back=True),
            }
            blockers, warnings, optional = [], [], []

            def block(field, message, action):
                blockers.append({"field": field, "message": message, "action": action})

            def warn(field, message, action):
                warnings.append({"field": field, "message": message, "action": action})

            configured_roles = {
                "condition_col": [audit_raw.get("condition_col"), deg_raw.get("condition_col"), abundance_raw.get("condition_col")],
                "donor_col": [audit_raw.get("donor_col"), deg_raw.get("donor_col"), abundance_raw.get("donor_col")],
                "batch_col": [audit_raw.get("batch_col"), integration_raw.get("batch_col")],
            }
            effective_roles = {
                "condition_col": [
                    self._effective_column(values, "condition_col", "condition")
                    for capability in ("dataset_audit", "deg", "differential_abundance")
                    for values in all_parameters(capability)
                ],
                "donor_col": [
                    self._effective_column(values, "donor_col", fallback)
                    for capability, fallback in (("dataset_audit", audit_donor_fallback),
                                                 ("deg", stats_donor_fallback),
                                                 ("differential_abundance", stats_donor_fallback))
                    for values in all_parameters(capability)
                ],
                "batch_col": [
                    self._effective_column(values, "batch_col", "batch" if "batch" in names else "", empty_falls_back=True)
                    for capability in ("dataset_audit", "integration")
                    for values in all_parameters(capability)
                ],
            }
            for field, columns in configured_roles.items():
                configured = sorted({str(column).strip() for column in columns if column not in (None, "")})
                for column in configured:
                    if column not in names:
                        block(field, f"最终计划引用的 {field} 列「{column}」不在当前文件中。", "从当前文件的列下拉框重新选择，或修改高级 JSON 的对应字段。")
                actual = sorted({column for column in effective_roles[field] if column})
                if len(actual) > 1:
                    block(field, "审计与分析任务的最终列映射不一致：" + "、".join(actual), "统一表单映射或高级 JSON 中的该字段；高级 JSON 覆盖表单。")

            if not any(configured_roles["condition_col"]) and "condition" in names:
                block("condition_col", "文件含有 condition 列，但尚未确认它是否代表实验条件。", "在实验设计中明确选择条件列；专家 JSON 中的列映射也会纳入检查。")
            donor_candidates = sorted({audit_donor_fallback, stats_donor_fallback} - {""})
            if not any(configured_roles["donor_col"]) and donor_candidates:
                block("donor_col", f"文件含有「{'、'.join(donor_candidates)}」候选列；后端方法默认顺序可能选用不同列。", "确认它是否代表独立生物学重复后，从供体列下拉框明确映射。")
            if not any(configured_roles["batch_col"]) and "batch" in names:
                block("batch_col", "文件含有 batch 列，但尚未确认它代表需要使用的批次。", "确认批次含义后明确映射，或在专家 JSON 中明确设置。")
            if form.sample_col and form.sample_col not in names:
                block("sample_col", f"样本标识列「{form.sample_col}」不在当前文件中。", "从当前文件的列下拉框重新选择。")

            condition_col = effective["audit_condition"]
            summary_donor = effective["deg_donor"] or effective["audit_donor"]
            summary_batch = effective["integration_batch"] or effective["audit_batch"]
            condition_values, condition_missing, condition_groups = [], None, []
            if condition_col in names:
                series = obs[condition_col]
                condition_missing = int(series.isna().sum())
                condition_values = [str(value) for value in series.dropna().unique().tolist()]
                for value in condition_values[:100]:
                    mask = series.astype(str).eq(value) & series.notna()
                    group = {"condition": value, "cells": int(mask.sum())}
                    if form.sample_col and form.sample_col in names:
                        group["samples"] = int(obs.loc[mask, form.sample_col].nunique(dropna=True))
                    if summary_donor and summary_donor in names:
                        group["donors"] = int(obs.loc[mask, summary_donor].nunique(dropna=True))
                    condition_groups.append(group)

            purpose = form.research_purpose
            deg_a = deg.get("condition_a") if plan else deg.get("condition_a", form.condition_a)
            deg_b = deg.get("condition_b") if plan else deg.get("condition_b", form.condition_b)
            abundance_a = abundance.get("condition_a") if plan else abundance.get("condition_a", form.condition_a)
            abundance_b = abundance.get("condition_b") if plan else abundance.get("condition_b", form.condition_b)
            primary_a, primary_b = (abundance_a, abundance_b) if purpose == "abundance" else (deg_a, deg_b)
            selected_a = str(primary_a) if primary_a not in (None, "") else ""
            selected_b = str(primary_b) if primary_b not in (None, "") else ""
            deg_a = str(deg_a) if deg_a not in (None, "") else ""
            deg_b = str(deg_b) if deg_b not in (None, "") else ""
            abundance_a = str(abundance_a) if abundance_a not in (None, "") else ""
            abundance_b = str(abundance_b) if abundance_b not in (None, "") else ""
            if bool(selected_a) != bool(selected_b):
                block("condition_a", "最终计划中的条件 A 和 B 必须同时填写。", "补齐另一个条件值，或从专家 JSON 中移除单独一侧。")
            if selected_a and selected_b and selected_a == selected_b:
                block("condition_b", "最终计划中的条件 A 和 B 不能相同。", "选择两个不同的实际条件值。")
            for capability, task_params, col_key in (
                ("deg", deg, "deg_condition"), ("differential_abundance", abundance, "abundance_condition")
            ):
                a = task_params.get("condition_a") if plan else task_params.get("condition_a", form.condition_a)
                b = task_params.get("condition_b") if plan else task_params.get("condition_b", form.condition_b)
                col = effective[col_key]
                if bool(a) != bool(b):
                    block("condition_a", f"{capability} 的最终任务只配置了一个条件值。", "为该任务补齐 A/B 条件，或移除专家 JSON 中不完整的对比。")
                if a not in (None, "") and b not in (None, ""):
                    if col not in names:
                        block("condition_col", f"{capability} 的条件列「{col}」不可用。", "选择当前文件中实际存在的条件列。")
                    else:
                        actual_values = {str(value) for value in obs[col].dropna().unique()}
                        for value in (a, b):
                            if str(value) not in actual_values:
                                block("condition_a" if str(value) == str(a) else "condition_b",
                                      f"条件值「{value}」不在最终映射列「{col}」中。", "从条件值下拉框选择当前文件中实际存在的值，或更正专家 JSON。")

            if purpose in {"deg", "abundance"} and not (selected_a and selected_b):
                block("condition_a", "此研究问题需要明确指定两个不同的条件值。", "选择条件 A 和 B；如果由专家 JSON 设置，请检查最终计划。")
            if condition_col in names and len(condition_values) == 2 and not (selected_a and selected_b):
                block("condition_a", "后端会在恰好两个条件值时自动建立比较；请先确认比较方向。", "在条件 A/B 下拉框明确选择两组，避免依赖自动顺序。")

            stats_donor = effective["deg_donor"]
            sample_count = int(obs[form.sample_col].nunique(dropna=True)) if form.sample_col in names and form.sample_col else None
            summary_donor_count = int(obs[summary_donor].nunique(dropna=True)) if summary_donor in names and summary_donor else None
            donor_by_condition = {row["condition"]: row.get("donors") for row in condition_groups}
            donors_a, donors_b = set(), set()
            if stats_donor in names and effective["deg_condition"] in names and deg_a and deg_b:
                conds = obs[effective["deg_condition"]].astype(str)
                donor_values = obs[stats_donor].astype("string")
                donors_a = set(donor_values.loc[conds.eq(deg_a)].dropna().tolist())
                donors_b = set(donor_values.loc[conds.eq(deg_b)].dropna().tolist())
            overlap = donors_a & donors_b
            abundance_donors_a, abundance_donors_b = set(), set()
            abundance_donor = effective["abundance_donor"]
            abundance_condition = effective["abundance_condition"]
            if abundance_donor in names and abundance_condition in names and abundance_a and abundance_b:
                abundance_conds = obs[abundance_condition].astype(str)
                abundance_values = obs[abundance_donor].astype("string")
                abundance_donors_a = set(abundance_values.loc[abundance_conds.eq(abundance_a)].dropna().tolist())
                abundance_donors_b = set(abundance_values.loc[abundance_conds.eq(abundance_b)].dropna().tolist())
            abundance_overlap = abundance_donors_a & abundance_donors_b
            method = tasks.get("deg", {}).get("method") or config.get("method_overrides", {}).get("deg")
            is_pydeseq = str(method or "").startswith("pydeseq2") or (not plan and bool(config.get("advanced_analysis")))
            paired = bool(deg.get("paired", deg.get("paired_design", False)))
            try:
                min_donors = max(2, int(deg.get("min_donors", 2)))
            except (TypeError, ValueError):
                min_donors = 2
                block("min_donors", "最终 PyDESeq2 最低供体数不是有效整数。", "设置至少 2 个供体，或修正专家 JSON。")
            if deg.get("min_donors") not in (None, ""):
                try:
                    if int(deg["min_donors"]) < 2:
                        block("min_donors", "最终 PyDESeq2 最低供体数不能少于 2。", "设置至少 2 个供体。")
                except (TypeError, ValueError):
                    pass
            if is_pydeseq:
                donor_series = obs[stats_donor] if stats_donor in names and stats_donor else None
                donor_complete = bool(donor_series is not None and not donor_series.isna().any()
                                      and donor_series.astype(str).str.strip().ne("").all())
                if not donor_complete:
                    block("donor_col", "PyDESeq2 需要完整的供体 ID；当前映射缺失、含空值或无法解析。", "选择完整的供体列，并检查每个供体的条件观测。")
                if paired and (not deg_a or not deg_b or len(overlap) < min_donors or donors_a != donors_b):
                    block("paired", f"配对设计要求至少 {min_donors} 个供体都同时出现在 A、B 两组；当前共同供体为 {len(overlap)} 个。", "检查供体与条件映射，确保同一供体在两个条件均有观测且无缺失。")
                elif not paired and overlap:
                    block("paired", f"独立组 PyDESeq2 中有 {len(overlap)} 个供体同时出现在 A、B 两组。", "若实验为同一供体配对设计，启用配对；否则修正供体或条件映射。")
                if deg_a and deg_b and min(len(donors_a), len(donors_b)) < min_donors:
                    block("min_donors", f"PyDESeq2 要求每个条件至少 {min_donors} 个供体；当前为 {len(donors_a)} 与 {len(donors_b)}。", "核对供体映射；细胞数不能替代生物学重复数。")
                counts_layer = deg.get("counts_layer", "counts")
                if counts_layer != "counts" and counts_layer not in data.layers:
                    block("counts_layer", f"PyDESeq2 指定的计数层「{counts_layer}」不存在。", "选择当前 AnnData 中存在的原始整数 counts layer。")
                if counts_layer == "counts" and counts_layer not in data.layers:
                    warn("counts_layer", "源文件没有 counts layer；标准预处理可能从 X 建立 counts layer，矩阵内容仍需执行时验证。", "确认 X 的含义；运行阶段会检查原始计数是否有效。")
                optional.append({"field": "counts_layer", "message": "预检没有逐项读取矩阵以验证非负整数计数；PyDESeq2 会在执行时检查。"})
            elif deg_a and deg_b and overlap and not paired:
                block("paired", f"独立组差异表达中有 {len(overlap)} 个供体同时出现在 A、B 两组。", "若实验为同一供体配对设计，启用高级统计和配对；否则修正供体或条件映射。")
            elif deg_a and deg_b and min(len(donors_a), len(donors_b)) < 2:
                warn("donor_col", "每组尚未显示至少两个供体；差异表达可能退化为单细胞探索性结果。", "检查供体映射与样本设计；不要把细胞数量当作生物学重复。")
            if paired and not is_pydeseq:
                block("paired", "配对设计只支持已注册的 PyDESeq2 差异表达方法。", "启用高级统计并在最终计划中选择 PyDESeq2，或取消配对选项。")

            if purpose == "abundance":
                if not abundance_donor or abundance_donor not in names:
                    block("donor_col", "细胞丰度比较需要明确的供体/生物学重复列。", "从供体列下拉框映射独立重复单位。")
                if min(len(abundance_donors_a), len(abundance_donors_b)) < 2:
                    block("donor_col", "当前 A/B 映射未显示每组至少两个供体；数据审计可能省略丰度比较。", "确认每个供体至少有一个样本，并检查条件与供体列映射。")
                if abundance_overlap:
                    block("paired", "当前供体同时出现在两个条件；现有丰度方法只支持独立供体。", "为丰度比较选择独立设计；配对丰度检验尚未由后端支持。")
            if purpose == "trajectory":
                root_id = trajectory.get("root_cell_id") if plan else trajectory.get("root_cell_id", form.root_cell_id)
                if not root_id:
                    block("root_cell_id", "拟时序研究问题需要一个有依据的根细胞 ID。", "填写当前文件中具有生物学依据的细胞 ID；不要随机选择。")
                elif str(root_id) not in {str(value) for value in obs.index}:
                    block("root_cell_id", f"根细胞 ID「{root_id}」不在当前文件的 obs 索引中。", "从原始数据的细胞 ID 中选择有效根细胞。")
            if purpose == "batch" and not effective["integration_batch"]:
                block("batch_col", "批次研究问题需要明确的批次列。", "从元数据列下拉框选择批次列。")
            extension = (config.get("analysis_extensions") or {}).get("cell_annotation")
            annotation_params = parameters("cell_annotation")
            if not tasks_by_capability.get("cell_annotation"):
                annotation_params = dict(extension) if isinstance(extension, dict) else {}
                annotation_params.update(raw_parameters.get("cell_annotation") or {})
            if purpose == "annotation":
                if plan and not tasks_by_capability.get("cell_annotation"):
                    block("celltypist_model_path", "最终计划没有包含参考模型注释任务。", "检查专家 JSON 中的 analysis_extensions.cell_annotation；该任务必须未被关闭。")
                model_path = annotation_params.get("model_path")
                if not model_path or not Path(model_path).is_file():
                    block("celltypist_model_path", "本地 CellTypist 注释需要一个存在的模型文件。", "填写可读取的本地模型路径；WebUI 不会下载模型。")
            if annotation_params.get("use_as_cell_type") and not annotation_params.get("model_path"):
                block("celltypist_model_path", "选择把预测标签用于后续分析时，必须同时配置本地模型文件。", "填写模型路径，或取消此选项。")

            for field, actual in (("condition_col", effective["audit_condition"]), ("donor_col", summary_donor),
                                  ("batch_col", effective["audit_batch"]), ("condition_a", selected_a),
                                  ("condition_b", selected_b), ("paired", paired)):
                entered = getattr(form, field)
                if str(entered).strip() != str(actual).strip() and (entered or actual):
                    warn(field, f"最终计划使用「{actual or '未设置'}」，与表单值不同。", "高级 JSON 优先；检查最终配置。")
            optional.append({"field": "data", "message": "此预检检查 AnnData 结构与元数据；不会验证表达矩阵质量或生物学合理性。"})
            multiple_targets = len(form.target_cell_types) > 1
            if multiple_targets:
                scope = "重复数摘要按全体细胞计算，未逐个目标细胞分支核对细胞数与供体覆盖；这里的预检通过不代表每个目标分支均满足运行条件。运行时数据审计仍会调整或省略不满足条件的分支。"
                optional.append({"field": "target_cell_types", "message": scope})
            else:
                scope = "此预检检查最终计划的元数据级输入；矩阵内容与运行时数据审计仍可能调整或省略不满足条件的分支。"
            overview = {
                "cells": int(data.n_obs), "genes": int(data.n_vars),
                "samples": {"column": form.sample_col or None, "count": sample_count,
                            "message": "样本标识唯一值仅作样本概览，不代表独立生物学重复。"},
                "donors": {"column": summary_donor or None, "count": summary_donor_count,
                           "confirmed": summary_donor in {str(value).strip() for value in configured_roles["donor_col"] if value not in (None, "")},
                           "by_condition": donor_by_condition, "paired_overlap": len(overlap)},
                "conditions": {"column": condition_col or None, "count": len(condition_values),
                               "confirmed": condition_col in {str(value).strip() for value in configured_roles["condition_col"] if value not in (None, "")},
                               "selected_a": selected_a or None, "selected_b": selected_b or None,
                               "missing": condition_missing, "groups": condition_groups,
                               "truncated": len(condition_values) > 100},
                "batches": {"column": summary_batch or None,
                            "confirmed": summary_batch in {str(value).strip() for value in configured_roles["batch_col"] if value not in (None, "")},
                            "count": int(obs[summary_batch].nunique(dropna=True)) if summary_batch in names and summary_batch else None,
                            "missing": int(obs[summary_batch].isna().sum()) if summary_batch in names and summary_batch else None},
                "scope": scope,
            }
            return {"overview": overview,
                    "checklist": {"blockers": blockers, "warnings": warnings, "optional": optional},
                    "can_run": not blockers}
        finally:
            data.file.close()

    def preview(self, payload):
        from eacbp.orchestrator.planning import preview_study_plan
        form, data, manifest, config = self.prepare(payload)
        plan = preview_study_plan(manifest, config)
        design = self.design({**form.model_dump(), "data": str(data)}, plan=plan, prepared_config=config)
        manifest.experimental_design.total_samples = design["overview"]["samples"]["count"] or 0
        required = {"anndata"}
        by_method = {"harmonypy_v1": ["harmonypy", "scanpy"],
                     "scanpy_leiden_umap_v1": ["scanpy", "leidenalg"],
                     "scanpy_dpt_v1": ["scanpy"], "cellrank_fate_v1": ["cellrank"],
                     "pydeseq2_pseudobulk_v1": ["pydeseq2"],
                     "pydeseq2_leave_one_donor_out_v1": ["pydeseq2"],
                     "decoupler_ulm_v2": ["decoupler"],
                     "celltypist_local_v1": ["celltypist"]}
        for task in plan["tasks"]:
            required.update(by_method.get(task.get("method"), []))
        installed = self.settings()["packages"]
        missing = sorted(package for package in required if not installed.get(package))
        return {"form": form.model_dump(), "data": str(data),
                "manifest": manifest.model_dump(mode="json"), "config": config,
                "plan": plan, "missing_packages": missing, "design": design,
                "can_run": plan["valid"] and not missing and design["can_run"]}

    def start(self, payload):
        with self.directory_guard:
            return self._start(payload)

    def _start(self, payload):
        prepared = self.preview(payload)
        if not prepared["can_run"]:
            raise ValueError("计划包含错误或缺少依赖，请先修复预览中的问题")
        job = self.manager.submit(
            "run", run_dir=self.manager.new_run_dir(prepared["manifest"]["study_id"]),
            manifest=prepared["manifest"], config=prepared["config"], data=prepared["data"],
            webui_context={"research_purpose": prepared["form"]["research_purpose"]},
        )
        return {"job": job, "run_id": Path(job["run_dir"]).relative_to(self.manager.runs_dir).as_posix()}

    def results(self, run_id):
        from .results import empty_results, result_views

        run = self.manager.run_path(run_id)
        if not run.is_dir():
            raise FileNotFoundError("运行记录不存在")
        if any(job.get("run_dir") == str(run) and job["status"] in {"queued", "running"}
               for job in self.manager.jobs()):
            return empty_results("后台操作尚未完成，请完成后刷新结果预览。")
        return result_views(run)

    def deg_table(self, run_id, *, task_id, page=1, size=50, query="",
                  fdr_max=None, significant_only=False, abs_log2fc_min=None):
        from .results import deg_table_page

        run = self.manager.run_path(run_id)
        if not run.is_dir():
            raise FileNotFoundError("运行记录不存在")
        latest = next((job for job in self.manager.jobs() if job.get("run_dir") == str(run)), None)
        if latest and latest["status"] in {"queued", "running", "unknown"}:
            raise BlockingIOError("该运行仍在执行或状态待检查，完成后再查询结果表。")
        if not task_id or len(task_id) > 200:
            raise ValueError("请提供有效的差异表达任务编号。")
        return deg_table_page(run, task_id, page=page, size=size, query=query,
                              fdr_max=fdr_max, significant_only=significant_only,
                              abs_log2fc_min=abs_log2fc_min)

    def deg_tables(self, run_id):
        from .results import deg_table_catalog

        run = self.manager.run_path(run_id)
        if not run.is_dir():
            raise FileNotFoundError("运行记录不存在")
        latest = next((job for job in self.manager.jobs() if job.get("run_dir") == str(run)), None)
        if latest and latest["status"] in {"queued", "running", "unknown"}:
            raise BlockingIOError("该运行仍在执行或状态待检查，完成后再读取结果目录。")
        return deg_table_catalog(run)

    def deg_table_csv(self, run_id, *, task_id):
        from .results import audited_deg_csv

        run = self.manager.run_path(run_id)
        if not run.is_dir():
            raise FileNotFoundError("运行记录不存在")
        latest = next((job for job in self.manager.jobs() if job.get("run_dir") == str(run)), None)
        if latest and latest["status"] in {"queued", "running", "unknown"}:
            raise BlockingIOError("该运行仍在执行或状态待检查，完成后再导出结果表。")
        if not task_id or len(task_id) > 200:
            raise ValueError("请提供有效的差异表达任务编号。")
        return audited_deg_csv(run, task_id)

    def export_bundle(self, run_id):
        """Export a completed, currently verified study bundle outside its run tree."""
        from .results import _verified_registry
        from eacbp.application.study_service import export_study

        run = self.manager.run_path(run_id)
        if not run.is_dir():
            raise FileNotFoundError("运行记录不存在")
        latest = next((job for job in self.manager.jobs() if job.get("run_dir") == str(run)), None)
        if latest and latest["status"] in {"queued", "running", "unknown"}:
            raise BlockingIOError("该运行仍在执行或状态待检查，完成后再导出证据包。")
        saved = read_json(run / "run_config.json", {})
        if saved.get("status") != "success" or saved.get("import_completed") is not True:
            raise ValueError("只有已成功完成数据导入的运行才能导出证据包。")
        snapshot, registry, _artifact_root = _verified_registry(run)
        contexts_by_task = {}
        for context in snapshot.audit_contexts.values():
            contexts_by_task.setdefault(context.task_id, []).append(context)
        for task in snapshot.task_results:
            contexts = contexts_by_task.get(task.task_id, [])
            if task.status.value != "success":
                raise ValueError(f"快照任务 {task.task_id} 状态为 {task.status.value}；完整证据包仅接受全部成功的运行。")
            if len(contexts) != 1:
                raise ValueError(f"任务 {task.task_id} 缺少唯一审计上下文，不能导出为完整有效证据包。")
            context = contexts[0]
            record = registry.get_audit_record(context.signature, context.receipt.get("contract")) if context.receipt else None
            if (not context.receipt_present or context.status != "passed" or not context.overall_passed
                    or context.stop_rule_triggered or record is None or record.status.value != "passed"
                    or not record.overall_passed or record.stop_rule_triggered
                    or record.task_id != task.task_id
                    or record.auditor_version != context.auditor_version
                    or record.auditor_fingerprint != context.auditor_fingerprint):
                raise ValueError(f"任务 {task.task_id} 没有与快照一致的通过审计回执，不能导出为完整有效证据包。")
            if not context.output_artifacts:
                raise ValueError(f"任务 {task.task_id} 没有登记输出产物，不能导出为完整有效证据包。")
            registry.verify_many(context.output_artifacts)
        export_root = contained(self.manager.workspace, self.manager.workspace / ".eacbp" / "webui_exports")
        export_root.mkdir(parents=True, exist_ok=True)
        output = export_root / f"{run.name}-{uuid4().hex[:10]}.zip"
        result = export_study(run_dir=run, output=output)
        if Path(result["bundle"]).resolve() != output.resolve() or not output.is_file():
            raise ValueError("证据包服务没有生成预期的受控下载文件。")
        return output, output.name

    def methods_report(self, run_id):
        """Human-readable provenance from the saved, currently verified snapshot."""
        from .results import _verified_registry

        run = self.manager.run_path(run_id)
        if not run.is_dir():
            raise FileNotFoundError("运行记录不存在")
        latest = next((job for job in self.manager.jobs() if job.get("run_dir") == str(run)), None)
        if latest and latest["status"] in {"queued", "running", "unknown"}:
            raise BlockingIOError("该运行仍在执行或状态待检查，完成后再下载方法与参数记录。")
        snapshot, registry, _artifact_root = _verified_registry(run)
        contexts = {}
        for context in snapshot.audit_contexts.values():
            contexts.setdefault(context.task_id, []).append(context)
        lines = [f"# 方法、参数与软件版本：{snapshot.manifest.study_id}", "",
                 f"历史研究：{snapshot.manifest.title}", "",
                 "本文件由已保存的证据快照生成。它记录历史运行，不重算分析、不重新审计科学结果；每项只列出快照中成功且有唯一通过审计回执的任务。完整结果与证据产物请下载完整证据包。", ""]
        simulated = (snapshot.config.get("mode") != "real"
                     or any(item.summary_metrics.get("is_simulated") for item in snapshot.artifact_metadata.values())
                     or any(item.is_simulated for item in snapshot.evidence_graph.evidence_nodes.values()))
        if simulated:
            lines.extend(["> 演示 / 合成证据标记：此运行包含模拟数据或模拟证据，不代表真实生物学观测。", ""])
        admitted = []
        for task in snapshot.task_results:
            matches = contexts.get(task.task_id, [])
            if task.status.value != "success" or len(matches) != 1:
                continue
            context = matches[0]
            record = registry.get_audit_record(context.signature, context.receipt.get("contract")) if context.receipt else None
            if (not context.receipt_present or context.status != "passed" or not context.overall_passed
                    or context.stop_rule_triggered or record is None or record.status.value != "passed"
                    or not record.overall_passed or record.stop_rule_triggered or record.task_id != task.task_id
                    or record.auditor_version != context.auditor_version
                    or record.auditor_fingerprint != context.auditor_fingerprint):
                continue
            metadata = registry.get_metadata_many(context.output_artifacts)
            registry.verify_many(context.output_artifacts)
            versions = {}
            artifacts = []
            for uri in context.output_artifacts:
                item = metadata[uri]
                versions.update(item.software_versions or {})
                artifacts.append({"uri": uri, "type": item.type.value,
                                  "parameters": item.parameters, "software_versions": item.software_versions})
            admitted.append((task, context, versions, artifacts))
        saved_parameters = snapshot.config.get("capability_parameters", {}) or {}
        for task, context, versions, artifacts in admitted:
            lines.extend([f"## {task.task_id} · {task.capability}", "",
                          f"- 状态：{task.status.value}",
                          f"- 方法：{task.method_used or '快照未记录'}",
                          f"- 审计：{context.status}；通过={context.overall_passed}；停止规则={context.stop_rule_triggered}",
                          "- 软件版本（历史产物元数据）：", ""])
            if versions:
                lines.extend([f"  - {name}: {value}" for name, value in sorted(versions.items())])
            else:
                lines.append("  - 快照产物未记录软件版本。")
            lines.extend(["", "任务参数与指标：", "", "```json",
                          json.dumps({"snapshot_capability_parameters": saved_parameters.get(task.capability, {}),
                                      "artifact_parameters": {item["uri"]: item["parameters"] for item in artifacts},
                                      "task_metrics": task.metrics}, ensure_ascii=False,
                                     indent=2, default=str), "```", "", "输出产物记录：", "", "```json",
                          json.dumps(artifacts, ensure_ascii=False, indent=2, default=str), "```", ""])
        if not admitted:
            lines.extend(["没有找到可列出的成功且通过审计的任务。", ""])
        return "\n".join(lines)

    def snapshot_download(self, run_id):
        """Return the historical snapshot only after its registry scope is revalidated."""
        from .results import _verified_registry, MAX_JSON_BYTES, _bounded_file

        run = self.manager.run_path(run_id)
        if not run.is_dir():
            raise FileNotFoundError("运行记录不存在")
        latest = next((job for job in self.manager.jobs() if job.get("run_dir") == str(run)), None)
        if latest and latest["status"] in {"queued", "running", "unknown"}:
            raise BlockingIOError("该运行仍在执行或状态待检查，完成后再下载方法与软件版本快照。")
        snapshot, _registry, _artifact_root = _verified_registry(run)
        path, _size = _bounded_file(run, run / "snapshot.json", MAX_JSON_BYTES)
        return path, f"{snapshot.manifest.study_id}-snapshot.json"
