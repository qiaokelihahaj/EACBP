"""Review pinned, externally computed descriptive tables without re-running scRNA.

The frozen_scoring_tables_v1 format describes small exports, not a native EACBP
study journal. Passing this review cannot admit native statistical evidence.
"""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path
import re
from typing import Any

import pandas as pd

from eacbp._atomic_files import publish_new_text
from eacbp.application.descriptive_checks import check_descriptive_tables


_ROLES = {"cells", "gene_summaries", "module_summaries", "scoring", "sensitivity", "modules"}
_MAX_BYTES = 64 * 1024 * 1024


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _pin(value: Any) -> str:
    _require(isinstance(value, str) and re.fullmatch(r"[0-9a-f]{64}", value) is not None,
             "Expected an explicit lowercase SHA256 pin")
    return value


def _read_pinned(path: Path, expected: str) -> bytes:
    _pin(expected)
    with path.open("rb") as stream:
        data = stream.read(_MAX_BYTES + 1)
    _require(len(data) <= _MAX_BYTES, f"Small-table review input exceeds 64 MiB: {path}")
    _require(hashlib.sha256(data).hexdigest() == expected, f"SHA256 mismatch: {path}")
    return data


def _json(data: bytes) -> dict:
    def reject(value):
        raise ValueError(f"Nonfinite JSON value: {value}")

    def unique_pairs(pairs):
        result = {}
        for key, value in pairs:
            _require(key not in result, f"Duplicate JSON key: {key}")
            result[key] = value
        return result

    result = json.loads(data.decode("utf-8-sig"), parse_constant=reject, object_pairs_hook=unique_pairs)
    _require(isinstance(result, dict), "Expected a JSON object")
    return result


def _table(data: bytes) -> pd.DataFrame:
    # Preserve identifiers such as state '0'; numerical blanks remain missing.
    import csv
    header = next(csv.reader(io.StringIO(data.decode("utf-8-sig"))), [])
    _require(bool(header) and len(header) == len(set(header)), "CSV headers must be nonempty and unique")
    identifiers = {"cell_id", "sample_id", "cluster", "group_scope", "group_label", "module_id",
                   "mask", "reference_sample", "query_sample", "gene_symbol", "epithelial_subtype"}
    return pd.read_csv(io.BytesIO(data), dtype={key: str for key in identifiers.intersection(header)})


def review_descriptive_bundle(bundle: str | Path, expected_bundle_sha256: str) -> dict:
    """Fail closed on source or table inconsistencies; perform no writes.

    Input paths are explicit local files relative to the bundle, or absolute.
    Only six declared data roles are read. Text is never used as instructions.
    """
    bundle = Path(bundle).resolve()
    spec = _json(_read_pinned(bundle, expected_bundle_sha256))
    _require(spec.get("schema_version") == 1 and spec.get("format") == "frozen_scoring_tables_v1",
             "Unsupported descriptive bundle schema/format")
    for key in ("study_id", "population", "reference_sample", "query_sample"):
        _require(isinstance(spec.get(key), str) and bool(spec[key].strip()), f"Missing nonempty {key}")
    _require(spec["reference_sample"] != spec["query_sample"], "Reference and query must differ")
    minimum = spec.get("min_cells")
    _require(type(minimum) is int and minimum >= 1, "min_cells must be a positive integer")
    design = spec.get("design")
    _require(isinstance(design, dict) and design.get("observation_unit") == "library"
             and design.get("biological_replication") in {"unverified", "unreplicated"},
             "This format requires explicit library units and unverified/unreplicated biological sampling")
    sources = spec.get("sources")
    _require(isinstance(sources, dict) and set(sources) == _ROLES,
             "Declare exactly cells, gene_summaries, module_summaries, scoring, sensitivity, modules")
    loaded, pins = {}, []
    for role in sorted(_ROLES):
        entry = sources[role]
        _require(isinstance(entry, dict) and isinstance(entry.get("path"), str) and entry["path"].strip(),
                 f"Invalid source path for {role}")
        path = (bundle.parent / entry["path"]).resolve()
        _require(path != bundle and all(str(path) != item["path"] for item in pins),
                 "Each source role must bind a distinct file, separate from the bundle")
        data = _read_pinned(path, entry.get("sha256"))
        loaded[role] = _table(data) if role in {"cells", "gene_summaries", "module_summaries"} else _json(data)
        pins.append({"role": role, "path": str(path), "sha256": entry["sha256"], "bytes": len(data)})

    scoring, sensitivity, definitions = (loaded[key] for key in ("scoring", "sensitivity", "modules"))
    source_sha = {item["role"]: item["sha256"] for item in pins}
    for name, manifest in (("scoring", scoring), ("sensitivity", sensitivity)):
        _require(manifest.get("status") == "complete" and manifest.get("inferential_statistics") is False,
                 f"{name} must declare complete, descriptive-only execution")
    _require(scoring.get("module_file_sha256") == source_sha["modules"]
             and sensitivity.get("source_scoring_module_sha256") == source_sha["modules"]
             and sensitivity.get("modules", {}).get("sha256") == source_sha["modules"],
             "Module definitions are not bound to the same scoring and sensitivity run")
    _require(sensitivity.get("source_scoring_manifest", {}).get("sha256") == source_sha["scoring"],
             "Sensitivity references a different scoring manifest")
    _require(scoring.get("module_set") == definitions.get("module_set_id")
             and scoring.get("module_version") == definitions.get("version"), "Module identity/version mismatch")
    parameters = sensitivity.get("parameters", {})
    for key in ("reference_sample", "query_sample", "min_cells"):
        _require(parameters.get(key) == spec[key], f"Sensitivity {key} differs from bundle")
    _require(parameters.get("state_column") == "cluster", "This export format requires the recorded cluster column")
    _require(sensitivity.get("original_scores_unchanged") is True, "Fixed-score sensitivity declaration missing")
    basis = sensitivity.get("source_score_basis", {})
    for key in ("normalization", "seed", "score_parameters", "score_background_gene_count"):
        _require(key in scoring and basis.get(key) == scoring[key], f"Changed/missing score basis: {key}")
    _require(scoring.get("n_epithelial_cells") == len(loaded["cells"])
             and sensitivity.get("n_source_cells") == len(loaded["cells"]), "Source cell totals disagree")
    _require(sensitivity.get("n_module_delta_rows") == len(loaded["module_summaries"]),
             "Module table has missing or extra rows relative to its execution manifest")
    module_records = definitions.get("modules")
    _require(isinstance(module_records, list) and all(isinstance(item, dict) and isinstance(item.get("id"), str)
                                                     for item in module_records), "Invalid module definitions")
    module_ids = [item["id"] for item in module_records]
    _require(module_ids and len(module_ids) == len(set(module_ids)), "Module IDs must be nonempty and unique")
    _require(set(loaded["module_summaries"]["module_id"]) == set(module_ids), "Module coverage differs from definitions")
    _require(set(sensitivity.get("score_columns", [])) == {"module_" + value for value in module_ids},
             "Score column declarations do not cover exactly the module definitions")

    checks = check_descriptive_tables(loaded["cells"], loaded["module_summaries"], loaded["gene_summaries"],
                                     reference_sample=spec["reference_sample"], query_sample=spec["query_sample"],
                                     population=spec["population"], min_cells=minimum)
    labels = {item["id"]: {key: item[key] for key in ("domain", "role", "label_zh") if key in item}
              for item in module_records}
    result = {
        "schema_version": 1, "review_status": "passed", "reviewed_utc": datetime.now(timezone.utc).isoformat(),
        "study_id": spec["study_id"], "population": spec["population"],
        "reference_sample": spec["reference_sample"], "query_sample": spec["query_sample"],
        "direction": "query_minus_reference", "min_cells": minimum, "design": design,
        "source_integrity": "passed", "raw_recomputation": "selected_exported_genes_only",
        "external_summary_consistency": "passed", "inference_status": "not_estimable",
        "native_evidence_admitted": False,
        "bundle": {"path": str(bundle), "sha256": expected_bundle_sha256}, "sources": pins,
        "reviewer_source_sha256": {
            path.name: hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (Path(__file__), Path(__file__).with_name("descriptive_checks.py"))},
        "module_definitions": labels, "checks": checks,
        "verification_scope": {
            "recomputed": "Selected gene UMI sums, detections and CPM using exported full-library total_counts; exact mask memberships",
            "checked": "Module summary arithmetic, group sizes and eligibility; hash-bound shared scoring declarations",
            "not_recomputed": ["FASTQ quantification", "full expression matrix", "gene mapping", "cell identities",
                               "module scores or control-gene background", "identity inclusion sensitivity", "gene ablation"],
            "independence": "All tables, states and masks derive from the same two libraries; none adds biological replication",
            "hash_scope": "Pins detect changes relative to supplied local records; they are not an external signature",
        },
    }
    # Guarantee a portable JSON receipt without implicit NaN/Infinity extensions.
    json.dumps(result, allow_nan=False)
    return result


def _cell(value: Any) -> str:
    return str(value).replace("|", "\\|").replace("\n", " ").replace("\r", " ").replace("<", "&lt;").replace(">", "&gt;")


def _number(value: Any) -> str:
    return "NA" if value is None else f"{value:.6g}" if isinstance(value, float) else _cell(value)


def format_descriptive_review(result: dict) -> str:
    """Always render the design boundary alongside numerical observations."""
    checks = result["checks"]
    counts = checks["checks"]
    primary_genes = [row for row in checks["raw_gene_summary"] if row["group_scope"] == "pooled" and row["mask"] == "all"]
    sample_counts = {row["sample_id"]: row["n_cells"] for row in primary_genes}
    lines = [f"# {_cell(result['study_id'])}：描述性结果复核", "",
             "来源校验和摘要一致性检查通过。下列差值方向为 query − reference。", "",
             f"reference：{_cell(result['reference_sample'])}；query：{_cell(result['query_sample'])}。",
             f"细胞数分别为 {sample_counts.get(result['reference_sample'], 'NA')} / {sample_counts.get(result['query_sample'], 'NA')}。",
             f"研究群体：{_cell(result['population'])}；观察单位为文库，独立生物学重复状态为 {_cell(result['design']['biological_replication'])}。",
             "本报告提供两文库之间的描述性观察，不估计处理效应，不输出组间 P 值、FDR 或置信区间。",
             "细胞、亚群及重叠过滤集合都属于同一来源，不能增加独立重复数。", "",
             "## 可复算范围", "",
             "逐细胞小表支持所导出基因的 UMI、检出率和 CPM 重算；CPM 使用全特征 total_counts 作为分母。",
             "模块检查覆盖摘要数值、样本量与同一评分背景的来源绑定；未重新计算模块分数或完整表达矩阵。",
             "上皮身份纳入敏感性与单基因消融仍未完成；过滤敏感性范围不代表置信区间。", ""]
    lines += ["## 核验结果", "",
              f"{counts['n_cells']} 个细胞、{counts['n_raw_genes']} 个导出基因：重算并核对 {counts['n_gene_rows_checked']} 行 UMI/检出汇总。",
              f"{counts['n_module_ids']} 个模块：核对 {counts['n_module_rows_checked']} 行摘要。",
              f"原基因汇总中其余 {counts['n_gene_rows_outside_selected_panel']} 行仅核验文件来源，未逐基因重算。",
              f"{counts['n_masks']} 个过滤规则对应 {counts['n_distinct_mask_cell_sets']} 个不同细胞集合；各集合仍然相互重叠。", "",
              "| 等价过滤规则 | reference 细胞数 | query 细胞数 |", "| --- | --- | --- |"]
    for row in checks["mask_equivalence"]:
        lines.append(f"| {_cell(', '.join(row['masks']))} | {row['n_cells_reference']} | {row['n_cells_query']} |")
    states = {(row["group_label"], row["sample_id"]): row["n_cells"] for row in checks["raw_gene_summary"]
              if row["group_scope"] == "subcluster_state" and row["mask"] == "all"}
    lines += ["", "## 状态内比较", "",
              f"以下为现有探索性表达状态。每库至少 {result['min_cells']} 个细胞是本项目的描述性展示门槛，不是独立重复要求或统计充分性证明。", "",
              "| 状态 | reference | query | 描述性比较资格 |", "| --- | --- | --- | --- |"]
    for state in sorted({key[0] for key in states}):
        a, b = (states.get((state, result[key]), 0) for key in ("reference_sample", "query_sample"))
        status = "缺失" if min(a, b) == 0 else "低细胞数" if min(a, b) < result["min_cells"] else "达到展示门槛"
        lines.append(f"| {_cell(state)} | {a} | {b} | {status} |")
    lines += ["", "合并群体的差异仍可能来自状态组成差异，不能据此声称同状态细胞自主效应。", ""]
    module_rows = [row for row in checks.get("module_sensitivity", []) if row["group_scope"] == "pooled"]
    if module_rows:
        lines += ["## 模块摘要", "", "这些数值描述原评分尺度内的差值，各模块的绝对值不用于排列机制强弱。", ""]
        lines += ["| 模块 | 全体差值 | 不同可比较集合数 | 最小差值 | 最大差值 | 方向 |", "| --- | --- | --- | --- | --- | --- |"]
        signs = {"consistent_positive": "均为正", "consistent_negative": "均为负", "consistent_zero": "均为零",
                 "sign_flip": "方向改变", "mixed_zero_and_nonzero": "含零值", "insufficient_distinct_comparable_cell_sets": "集合不足"}
        for row in module_rows:
            label = result["module_definitions"].get(row["module_id"], {}).get("label_zh", row["module_id"])
            lines.append(f"| {_cell(label)} | {_number(row['all_delta'])} | {row['n_distinct_comparable_cell_sets']} | {_number(row['delta_min'])} | {_number(row['delta_max'])} | {signs.get(row['sign_status'], _cell(row['sign_status']))} |")
    lines += ["", "## 逐细胞原始 UMI 重算", "", "检出率为至少一个 UMI 的细胞数 / 该库保留细胞数。这里只显示完整集合，全部过滤和状态结果见 JSON。", "",
              "| 基因 | reference 检出 | query 检出 | reference CPM | query CPM |", "| --- | --- | --- | --- | --- |"]
    for gene in sorted({row["gene_symbol"] for row in primary_genes}):
        rows = {row["sample_id"]: row for row in primary_genes if row["gene_symbol"] == gene}
        a, b = rows[result["reference_sample"]], rows[result["query_sample"]]
        lines.append(f"| {_cell(gene)} | {a['detected_cells']}/{a['n_cells']} | {b['detected_cells']}/{b['n_cells']} | {_number(a['pseudobulk_cpm'])} | {_number(b['pseudobulk_cpm'])} |")
    lines += ["", "## 来源与复现", "", f"证据包 SHA256：`{result['bundle']['sha256']}`。", "",
              "| 角色 | SHA256 |", "| --- | --- |"]
    lines += [f"| {_cell(item['role'])} | `{item['sha256']}` |" for item in result["sources"]]
    lines += ["", "完整路径、逐行重算和审核代码 SHA256 见同次 JSON 回执。",
              "这是外部描述性结果的限定范围复核，不等同于 EACBP 原生工作流重跑或机制证据准入。", "",
              "方法依据：[Confronting false discoveries in single-cell differential expression](https://www.nature.com/articles/s41467-021-25960-2)。", ""]
    return "\n".join(lines)


def publish_descriptive_review(result: dict, *, output: Path | None = None, report: Path | None = None) -> None:
    """Protect all sources and existing output versions before publishing."""
    destinations = [Path(path).resolve() for path in (output, report) if path is not None]
    protected = {Path(item["path"]).resolve() for item in [result["bundle"], *result["sources"]]}
    _require(len(destinations) == len(set(destinations)), "Review output paths must be distinct")
    _require(all(path not in protected and not path.exists() for path in destinations),
             "Review outputs must be new files, separate from source evidence")
    if output is not None:
        publish_new_text(output, json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n")
    if report is not None:
        publish_new_text(report, format_descriptive_review(result))
