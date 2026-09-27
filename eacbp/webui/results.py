"""Bounded, read-only chart views of outputs admitted in a saved snapshot."""
from __future__ import annotations

import math
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd

from .jobs import contained


MAX_FILE_BYTES = 32_000_000
MAX_TOTAL_BYTES = 128_000_000
MAX_HASH_BYTES = 2_000_000_000
MAX_JSON_BYTES = 8_000_000
MAX_POINTS = 2000
MAX_GENES = 200
MAX_DEG_TABLES = 12
MAX_DEG_TABLE_BYTES = 512_000_000
DEG_PAGE_SIZE = 50
MAX_DEG_PAGE_SIZE = 100
MAX_DEG_PAGE = 1_000_000
VERIFICATION_NOTICE = "来源文件哈希、快照上下文及已保存的通过审计回执已核对；未重新执行科学分析或审计。"


class PreviewSizeLimit(ValueError):
    """A valid artifact is too large for an in-memory browser projection."""


def empty_results(notice):
    return {"state": "unavailable", "notice": notice, "qc": None,
            "embedding": None, "deg": [], "warnings": [], "tasks": []}


def _bounded_file(root, path, limit=None):
    path = contained(root, path)
    if not path.is_file():
        raise ValueError("可视化来源文件不存在")
    size = path.stat().st_size
    if limit is not None and size > limit:
        raise PreviewSizeLimit("来源文件超过浏览器结果预览大小限制，请从运行目录查看")
    return path, size


def _check_expanded_size(path):
    """Reject large compressed matrices before the normal audited loader."""
    if path.suffix == ".npz":
        with ZipFile(path) as archive:
            total = sum(item.file_size for item in archive.infolist())
            if total > MAX_TOTAL_BYTES:
                raise PreviewSizeLimit("解压后的产物超过结果预览大小限制")
            return total
    elif path.suffix == ".h5ad":
        import h5py
        total = 0
        with h5py.File(path, "r") as handle:
            visited = set()

            def inspect(group):
                nonlocal total
                if group.id in visited:
                    return
                visited.add(group.id)
                for name in group:
                    link = group.get(name, getlink=True)
                    if not isinstance(link, h5py.HardLink):
                        raise ValueError("结果预览不读取 HDF5 外部或符号链接")
                    item = group[name]
                    if isinstance(item, h5py.Group):
                        inspect(item)
                    else:
                        if item.is_virtual or item.external:
                            raise ValueError("结果预览不读取 HDF5 外部数据")
                        total += int(item.size) * max(1, item.dtype.itemsize)
                        if total > MAX_TOTAL_BYTES:
                            raise PreviewSizeLimit("展开后的产物超过结果预览大小限制")
            inspect(handle)
        return total
    return path.stat().st_size


def _number(value):
    try:
        number = float(value)
    except (ValueError, TypeError):
        return None
    return number if math.isfinite(number) else None


def _label(value):
    return "未标注" if pd.isna(value) else str(value)[:160]


def _qc_view(metadata, payload, task):
    retained = getattr(payload, "n_obs", None)
    initial = _number(task.metrics.get("initial_cells"))
    # Counts derive from the admitted payload plus its saved task result.
    if retained is None or initial is None or initial < retained:
        return None
    return {"task_id": task.task_id, "source_uri": metadata.uri, "notice": VERIFICATION_NOTICE,
            "initial_cells": int(initial), "retained_cells": int(retained),
            "filtered_cells": int(initial - retained),
            "retention_rate": retained / initial if initial else 0}


def _embedding_view(metadata, payload, task):
    obsm = getattr(payload, "obsm", {})
    obs = getattr(payload, "obs", None)
    # Legacy baseline output may contain an X_umap alias. Its method is decisive.
    is_umap = task.method_used == "scanpy_leiden_umap_v1"
    key = "X_umap" if is_umap else "X_embedding_2d"
    if key not in obsm or obs is None:
        return None
    coords = np.asarray(obsm[key])
    if coords.ndim != 2 or coords.shape[1] < 2 or coords.shape[0] != len(obs):
        raise ValueError("二维坐标与细胞元数据不一致")
    finite = np.flatnonzero(np.isfinite(coords[:, :2]).all(axis=1))
    indexes = finite[np.linspace(0, len(finite) - 1, min(len(finite), MAX_POINTS), dtype=int)] if len(finite) else []
    clusters = obs["cluster"].map(_label) if "cluster" in obs else pd.Series("未标注", index=obs.index)
    types = obs["cell_type"].map(_label) if "cell_type" in obs else clusters
    counts = clusters.value_counts().sort_index()
    return {"task_id": task.task_id, "source_uri": metadata.uri, "notice": VERIFICATION_NOTICE,
            "method": task.method_used, "label": "UMAP" if is_umap else "二维展示投影（非 UMAP）",
            "total_cells": len(obs), "shown_cells": len(indexes), "sampled": len(indexes) < len(obs),
            "omitted_nonfinite": int(len(obs) - len(finite)),
            "points": [{"x": float(coords[i, 0]), "y": float(coords[i, 1]),
                        "cluster": str(clusters.iloc[i]), "cell_type": str(types.iloc[i])} for i in indexes],
            "groups": [{"label": str(label), "count": int(count)} for label, count in counts.items()]}


def _deg_view(metadata, payload, task):
    required = {"gene", "log2_fold_change", "fdr_q_value"}
    if not isinstance(payload, pd.DataFrame) or not required.issubset(payload.columns):
        return None
    table = payload.sort_values("fdr_q_value", kind="stable", na_position="last").head(MAX_GENES)
    context = {**metadata.parameters, **metadata.summary_metrics, **task.metrics}
    return {"task_id": task.task_id, "source_uri": metadata.uri, "notice": VERIFICATION_NOTICE,
            "method": task.method_used, "target_cell_type": context.get("target_cell_type"),
            "condition_a": context.get("condition_a"), "condition_b": context.get("condition_b"),
            "contrast_label": context.get("contrast_label"),
            "effect_definition": context.get("effect_definition"), "alpha": context.get("alpha", 0.05),
            "statistical_unit": context.get("statistical_unit", "unknown"),
            "total_genes": len(payload), "shown_genes": len(table),
            "rows": [{"gene": _label(row["gene"]),
                      "log2_fold_change": _number(row["log2_fold_change"]),
                      "fdr_q_value": _number(row["fdr_q_value"]),
                      "p_value": _number(row.get("p_value"))} for _, row in table.iterrows()]}


def result_views(run: Path):
    """Use the report's provenance gate; never infer admission from file names."""
    from eacbp.artifact.registry import ArtifactRegistry
    from eacbp.evidence.snapshot import load_study_snapshot, _verify_registry_scope

    response = empty_results("尚无完成的结果快照；分析结束后可查看结果。")
    if not (run / "snapshot.json").exists():
        return response
    try:
        snapshot_path, _ = _bounded_file(run, run / "snapshot.json", MAX_JSON_BYTES)
        artifact_root = contained(run, run / "artifacts")
        _bounded_file(run, artifact_root / ArtifactRegistry.INDEX_FILENAME, MAX_JSON_BYTES)
        contained(artifact_root, artifact_root / ".registry.lock")
        snapshot = load_study_snapshot(snapshot_path)
        registry = ArtifactRegistry(str(artifact_root))
        # Hash large ancestors in a stream without loading them into memory.
        total = 0
        for metadata in registry.list_artifacts():
            _, size = _bounded_file(artifact_root, metadata.storage_path)
            total += size
            if total > MAX_HASH_BYTES:
                raise PreviewSizeLimit("完整来源校验超过 2 GB 读取预算，请从运行目录查看")
        _verify_registry_scope(snapshot, registry)
        response["notice"] = VERIFICATION_NOTICE
        contexts = {}
        for context in snapshot.audit_contexts.values():
            contexts.setdefault(context.task_id, []).append(context)
        for task in snapshot.task_results:
            matching = contexts.get(task.task_id, [])
            audited = False
            audit_status = "missing_or_ambiguous"
            if len(matching) == 1:
                context = matching[0]
                audit_status = context.status
                try:
                    record = registry.get_audit_record(context.signature, context.receipt.get("contract")) if context.receipt else None
                    audited = bool(task.status.value == "success" and context.receipt_present
                                   and context.status == "passed" and context.overall_passed
                                   and not context.stop_rule_triggered and record is not None
                                   and record.task_id == task.task_id and record.status.value == "passed"
                                   and record.overall_passed and not record.stop_rule_triggered
                                   and record.auditor_version == context.auditor_version
                                   and record.auditor_fingerprint == context.auditor_fingerprint)
                except (KeyError, ValueError, RuntimeError, TypeError):
                    audited = False
            response["tasks"].append({"task_id": task.task_id, "capability": task.capability,
                                      "status": task.status.value, "method": task.method_used,
                                      "audit_status": audit_status, "audited": audited})
        for task in snapshot.task_results:
            if task.capability not in {"qc", "clustering", "deg"}:
                continue
            matching = contexts.get(task.task_id, [])
            # Ambiguous historical contexts are not silently promoted.
            if task.status.value != "success" or len(matching) != 1:
                response["warnings"].append(f"{task.task_id}：任务未成功或缺少唯一审计上下文，未展示。")
                continue
            context = matching[0]
            if not context.receipt_present or context.status != "passed" or not context.overall_passed or context.stop_rule_triggered:
                response["warnings"].append(f"{task.task_id}：缺少已通过的持久审计回执，未展示。")
                continue
            if task.capability == "deg" and len(response["deg"]) >= MAX_DEG_TABLES:
                response["warnings"].append("仅预览前 12 个差异表达结果表。")
                break
            try:
                expanded = 0
                for uri in context.output_artifacts:
                    metadata = registry.get_metadata(uri)
                    path, _ = _bounded_file(artifact_root, metadata.storage_path, MAX_FILE_BYTES)
                    expanded += _check_expanded_size(path)
                    if expanded > MAX_TOTAL_BYTES:
                        raise PreviewSizeLimit("任务输出展开后超过结果预览大小限制")
            except PreviewSizeLimit as exc:
                response["warnings"].append(f"{task.task_id}：{exc}；已跳过该任务预览。")
                continue
            values = registry.get_audited_many(context.signature, context.receipt.get("contract"),
                                               expected_auditor_fingerprint=context.auditor_fingerprint,
                                               expected_auditor_version=context.auditor_version)
            for metadata, payload in values:
                if task.capability == "qc":
                    response["qc"] = _qc_view(metadata, payload, task) or response["qc"]
                elif task.capability == "clustering":
                    response["embedding"] = _embedding_view(metadata, payload, task) or response["embedding"]
                elif task.capability == "deg":
                    view = _deg_view(metadata, payload, task)
                    if view:
                        response["deg"].append(view)
            del values
        if snapshot.config.get("mode") != "real" or any(
            metadata.summary_metrics.get("is_simulated") for metadata in snapshot.artifact_metadata.values()
        ) or any(node.is_simulated for node in snapshot.evidence_graph.evidence_nodes.values()):
            response["warnings"].insert(0, "演示 / 模拟数据，仅用于软件演示，不代表真实生物学观测。")
        if any(task.status.value != "success" for task in snapshot.task_results):
            response["warnings"].insert(0, "该研究存在未完成或未通过的任务；这里只展示已通过审计的局部结果。")
        response["state"] = "ready" if response["qc"] or response["embedding"] or response["deg"] else "unavailable"
        if response["state"] == "unavailable":
            response["notice"] = "当前快照没有可展示的 QC、聚类坐标或差异表达结果。"
        return response
    except (OSError, ValueError, RuntimeError, KeyError, TypeError, ImportError) as exc:
        # Discard all partial views if any live provenance check fails.
        return empty_results(f"结果预览不可用：{exc}")


def _verified_registry(run: Path):
    """Revalidate the complete snapshot scope before a table is queried/exported."""
    from eacbp.artifact.registry import ArtifactRegistry
    from eacbp.evidence.snapshot import load_study_snapshot, _verify_registry_scope

    snapshot_path, _ = _bounded_file(run, run / "snapshot.json", MAX_JSON_BYTES)
    artifact_root = contained(run, run / "artifacts")
    _bounded_file(run, artifact_root / ArtifactRegistry.INDEX_FILENAME, MAX_JSON_BYTES)
    contained(artifact_root, artifact_root / ".registry.lock")
    snapshot = load_study_snapshot(snapshot_path)
    registry = ArtifactRegistry(str(artifact_root))
    total = 0
    for metadata in registry.list_artifacts():
        _, size = _bounded_file(artifact_root, metadata.storage_path)
        total += size
        if total > MAX_HASH_BYTES:
            raise PreviewSizeLimit("完整来源校验超过 2 GB 读取预算，请从运行目录查看")
    _verify_registry_scope(snapshot, registry)
    return snapshot, registry, artifact_root


def _audited_deg_table(run: Path, task_id: str):
    """Return the one DEG TABLE source admitted by its current snapshot receipt."""
    snapshot, registry, artifact_root = _verified_registry(run)
    return _deg_table_metadata(snapshot, registry, artifact_root, task_id)


def _deg_table_metadata(snapshot, registry, artifact_root: Path, task_id: str):
    from eacbp.schemas.artifact import ArtifactType

    tasks = [task for task in snapshot.task_results if task.task_id == task_id]
    if len(tasks) != 1 or tasks[0].capability != "deg" or tasks[0].status.value != "success":
        raise ValueError("所选差异表达任务不存在或未成功；请从已完成结果中重新选择。")
    task = tasks[0]
    contexts = [context for context in snapshot.audit_contexts.values() if context.task_id == task_id]
    if len(contexts) != 1:
        raise ValueError("该差异表达任务没有唯一的快照审计上下文，不能查询或导出。")
    context = contexts[0]
    if (not context.receipt_present or context.status != "passed" or not context.overall_passed
            or context.stop_rule_triggered):
        raise ValueError("该差异表达任务没有通过审计的持久回执，不能查询或导出。")
    contract = context.receipt.get("contract")
    if not contract:
        raise ValueError("快照缺少该差异表达任务的审计契约，不能查询或导出。")
    record = registry.get_audit_record(context.signature, contract)
    if (record is None or record.task_id != task_id or record.status.value != "passed"
            or not record.overall_passed or record.stop_rule_triggered
            or record.auditor_version != context.auditor_version
            or record.auditor_fingerprint != context.auditor_fingerprint):
        raise ValueError("实时审计回执与快照不一致，不能查询或导出。")
    output_uris = list(context.output_artifacts)
    if not output_uris:
        raise ValueError("该差异表达任务没有登记输出产物。")
    # Stream-hash every sibling output. The table source is never trusted just
    # because its filename or task id looks plausible.
    metadata_by_uri = registry.get_metadata_many(output_uris)
    registry.verify_many(output_uris)
    candidates = [metadata_by_uri[uri] for uri in output_uris
                  if metadata_by_uri[uri].type == ArtifactType.TABLE]
    if len(candidates) != 1:
        raise ValueError("该差异表达任务没有唯一的完整结果表产物。")
    metadata = candidates[0]
    path, size = _bounded_file(artifact_root, metadata.storage_path, MAX_DEG_TABLE_BYTES)
    if size != metadata.size_bytes:
        raise ValueError("差异表达结果表大小与登记元数据不一致。")
    return task, context, metadata, path


def deg_table_catalog(run: Path):
    """List each full result table that has one live, passed audit context."""
    snapshot, registry, artifact_root = _verified_registry(run)
    catalog, warnings = [], []
    for task in snapshot.task_results:
        if task.capability != "deg" or task.status.value != "success":
            continue
        try:
            task, context, metadata, path = _deg_table_metadata(
                snapshot, registry, artifact_root, task.task_id)
            parameters = {**metadata.parameters, **metadata.summary_metrics, **task.metrics}
            catalog.append({"task_id": task.task_id, "source_uri": metadata.uri,
                            "method": task.method_used, "total_bytes": path.stat().st_size,
                            "target_cell_type": parameters.get("target_cell_type"),
                            "condition_a": parameters.get("condition_a"),
                            "condition_b": parameters.get("condition_b"),
                            "contrast_label": parameters.get("contrast_label"),
                            "effect_definition": parameters.get("effect_definition"), "alpha": parameters.get("alpha", 0.05),
                            "statistical_unit": parameters.get("statistical_unit", "unknown"),
                            "query_available": path.stat().st_size <= MAX_DEG_TABLE_BYTES,
                            "notice": VERIFICATION_NOTICE})
        except (OSError, ValueError, RuntimeError, KeyError, TypeError) as exc:
            warnings.append(f"{task.task_id}：未列为可查询结果表（{exc}）。")
    return {"tables": catalog, "warnings": warnings,
            "notice": "目录仅列出快照中成功且通过持久审计的差异表达 TABLE 产物；未重算分析。",
            "max_table_bytes": MAX_DEG_TABLE_BYTES}


def _table_chunk_rows(path: Path, *, page: int, size: int, query: str,
                      fdr_max: float | None, significant_only: bool,
                      abs_log2fc_min: float | None):
    import pandas as pd

    if not path.is_file():
        raise FileNotFoundError("差异表达结果表文件不存在。")
    result = []
    total = 0
    try:
        chunks = pd.read_csv(path, index_col=0, chunksize=4096, low_memory=False,
                             dtype={"gene": "string"}, keep_default_na=False)
        for chunk in chunks:
            chunk = chunk.reset_index(drop=True)
            required = {"gene", "log2_fold_change", "fdr_q_value"}
            if not required.issubset(chunk.columns):
                raise ValueError("结果表缺少 gene、log2_fold_change 或 fdr_q_value 列。")
            genes = chunk["gene"].fillna("").astype(str)
            mask = pd.Series(True, index=chunk.index)
            if query:
                mask &= genes.str.casefold().str.contains(query.casefold(), regex=False)
            fdr = pd.to_numeric(chunk["fdr_q_value"], errors="coerce")
            if significant_only:
                mask &= fdr.notna() & fdr.le(fdr_max)
            if abs_log2fc_min is not None:
                fold_change = pd.to_numeric(chunk["log2_fold_change"], errors="coerce").abs()
                mask &= fold_change.notna() & fold_change.ge(abs_log2fc_min)
            indexes = list(np.flatnonzero(mask.to_numpy()))
            first, last = (page - 1) * size, page * size
            overlap_start, overlap_end = max(first, total), min(last, total + len(indexes))
            if overlap_start < overlap_end:
                local_indexes = indexes[overlap_start - total:overlap_end - total]
                selected = chunk.iloc[local_indexes]
                for _, row in selected.iterrows():
                    def number(value):
                        try:
                            parsed = float(value)
                        except (TypeError, ValueError):
                            return None
                        return parsed if math.isfinite(parsed) else None
                    result.append({"gene": str(row["gene"]),
                                   "log2_fold_change": number(row.get("log2_fold_change")),
                                   "fdr_q_value": number(row.get("fdr_q_value")),
                                   "p_value": number(row.get("p_value"))})
            total += len(indexes)
    except (pd.errors.EmptyDataError, pd.errors.ParserError, UnicodeDecodeError) as exc:
        raise ValueError(f"无法读取差异表达 CSV 产物：{exc}") from exc
    return result, total


def deg_table_page(run: Path, task_id: str, *, page=1, size=DEG_PAGE_SIZE,
                   query="", fdr_max=None, significant_only=False, abs_log2fc_min=None):
    """Strictly paginated, searchable server-side view of an audited full DEG table."""
    if isinstance(page, bool) or not isinstance(page, int) or not 1 <= page <= MAX_DEG_PAGE:
        raise ValueError(f"page 必须是 1 到 {MAX_DEG_PAGE} 之间的整数。")
    if isinstance(size, bool) or not isinstance(size, int) or not 1 <= size <= MAX_DEG_PAGE_SIZE:
        raise ValueError(f"size 必须是 1 到 {MAX_DEG_PAGE_SIZE} 之间的整数。")
    if not isinstance(query, str) or len(query) > 200:
        raise ValueError("基因搜索词最多 200 个字符。")
    if not isinstance(significant_only, bool):
        raise ValueError("significant_only 必须是布尔值。")
    if fdr_max is not None:
        if isinstance(fdr_max, bool):
            raise ValueError("FDR 阈值必须是 0 到 1 之间的数字。")
        try:
            fdr_max = float(fdr_max)
        except (TypeError, ValueError) as exc:
            raise ValueError("FDR 阈值必须是 0 到 1 之间的数字。") from exc
        if not math.isfinite(fdr_max) or not 0 <= fdr_max <= 1:
            raise ValueError("FDR 阈值必须是 0 到 1 之间的数字。")
    if significant_only and fdr_max is None:
        raise ValueError("筛选显著基因时必须提供 FDR 阈值。")
    if abs_log2fc_min is not None:
        if isinstance(abs_log2fc_min, bool):
            raise ValueError("绝对 log2 fold change 阈值必须是非负有限数字。")
        try:
            abs_log2fc_min = float(abs_log2fc_min)
        except (TypeError, ValueError) as exc:
            raise ValueError("绝对 log2 fold change 阈值必须是非负有限数字。") from exc
        if not math.isfinite(abs_log2fc_min) or abs_log2fc_min < 0:
            raise ValueError("绝对 log2 fold change 阈值必须是非负有限数字。")
    task, context, metadata, path = _audited_deg_table(run, task_id)
    rows, total = _table_chunk_rows(path, page=page, size=size, query=query,
                                    fdr_max=fdr_max, significant_only=significant_only,
                                    abs_log2fc_min=abs_log2fc_min)
    parameters = {**metadata.parameters, **metadata.summary_metrics, **task.metrics}
    return {"task_id": task.task_id, "source_uri": metadata.uri,
            "method": task.method_used, "target_cell_type": parameters.get("target_cell_type"),
            "condition_a": parameters.get("condition_a"), "condition_b": parameters.get("condition_b"),
            "contrast_label": parameters.get("contrast_label"),
            "effect_definition": parameters.get("effect_definition"), "alpha": parameters.get("alpha", 0.05),
            "statistical_unit": parameters.get("statistical_unit", "unknown"),
            "notice": VERIFICATION_NOTICE, "page": page, "size": size,
            "total_items": total, "total_pages": math.ceil(total / size), "rows": rows}


def audited_deg_csv(run: Path, task_id: str):
    """Return a verified source path for streaming download of the full TABLE CSV."""
    _, _, metadata, path = _audited_deg_table(run, task_id)
    safe_task = "".join(char if char.isalnum() or char in "-_" else "_" for char in task_id)
    return path, f"deg-{safe_task}.csv"
