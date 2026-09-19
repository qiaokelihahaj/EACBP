# P3：统一运行入口与可复查结果

本轮落实 P3 的前两项：包内 CLI、只读计划预览，以及版本化的证据和报告快照。运行时限/GPU 策略执行、SCData 模块迁移和性能优化留待后续实施。

## 计划与执行共用合同解析

`eacbp.orchestrator.planning` 提供共享的 `build_study_tasks()` 和 `resolve_task_contract()`。编排器和预览都使用相同的扩展装配、依赖计算、方法路由、操作范围和能力声明校验，避免出现两套计划逻辑。重复产物 URI 在装配阶段拒绝。

`preview_study_plan(manifest, config)` 返回 JSON 兼容对象，列出任务、依赖、方法、输入输出、目标分支、缺失参数及静态错误。它不创建注册表或输出目录，不载入矩阵，也不运行任务计算。扩展规划器仍可能读取并哈希明确提供的资源文件。

预览发生在数据审计之前，结果标记为 `before_dataset_audit`：批次、重复数和条件信息可能导致实际运行选择不同方法或省略分析分支。DPT 缺少根细胞时提示其分支将被省略；FASTQ 定量缺少必要参考资源时标记阻断。预览不证明细胞 ID、输入内容、可选科学库或计算资源可用。

## 交付范围

- 安装后的统一入口包含计划、运行、恢复、日志检查、事务清理和报告重建。
- 运行目录保存实际使用的研究说明、配置和来源摘要，恢复使用同一运行目录与原有严格签名规则。
- 快照保存结构化证据、结论、任务结果、审计和来源元数据，不嵌入大矩阵。
- 报告重建需要验证来源完整性与当前审计状态；快照不能替代审计收据，也不能恢复已被撤销的科学通过状态。

## CLI 使用

安装带所需科学依赖的包后可使用 `eacbp`，源码环境也支持 `python -m eacbp`。原有 `scripts/` 入口保持兼容。新入口不依赖安装包以外的脚本目录。

```powershell
eacbp plan --manifest manifest.json --config config.json
eacbp run --manifest manifest.json --config config.json --data input.h5ad --run-dir outputs/run_001
eacbp resume --run-dir outputs/run_001
eacbp inspect --run-dir outputs/run_001
eacbp report --run-dir outputs/run_001 --output rebuilt_report.md
eacbp cleanup outputs/run_001/artifacts --older-than-days 30
```

`manifest.json` 使用 `StudyManifest` 字段，例如：

```json
{
  "study_id": "brain_study",
  "biological_design": {
    "species": "human",
    "tissue": "brain",
    "target_cell_types": ["Microglia", "Neurons"]
  }
}
```

`config.json` 使用 `RunConfig` 字段；可省略，实际运行按既有规则选择 profile。标准 DPT 的多目标根细胞等设置见 P2 文档。`plan` 的 JSON `valid=false` 或命令失败返回非零退出码。`cleanup` 默认只预览，删除需另加 `--apply`；阻断和错误同样返回非零。

新运行目录必须不存在，包括已有空目录也会拒绝。CLI 在导入前后核对原始 h5ad 哈希，将实际 manifest、规范化配置和来源摘要保存在 `run_config.json` 及辅助 JSON 文件。恢复使用注册表内已导入的副本，不重新读取或导入外部 h5ad；外部原文件后来移动或改变，不会悄悄改变该次研究。

CLI 使用运行目录锁覆盖计算、摘要、快照和报告发布，串行化同目录内的 CLI 操作。任务层原有 study 锁和严格恢复签名仍生效。源码/依赖或任务配置改变导致的严格恢复拒绝不会被 CLI 绕过。

每次正常返回的运行保存 `summaries/<run_id>.json` 和 `snapshots/<run_id>.json`，`summary.json` / `snapshot.json` 指向最新内容。交付快照或报告失败会使 CLI 返回失败，保留已计算产物与恢复日志。

## 快照校验边界

Python API 位于 `eacbp.evidence.snapshot`：

- `write_study_snapshot(...)` 检查任务、证据、审计、产物之间的关联及当前来源，再原子写入版本化 JSON。
- `load_study_snapshot(path)` 校验格式版本、各部分摘要、图引用和审计上下文，返回结构化 `StudySnapshot`。可通过 `rebuild_evidence_graph()` 获得独立证据图。
- `render_snapshot_report(path, artifact_registry)` 重新验证文件哈希、元数据、计算收据和当前审计状态后生成 Markdown，不重新运行科学计算。

所有接纳的计算证据都必须对应成功任务及持久化通过审计；仅有一份内存报告不够。原始导入可作为祖先保留，无需伪造审计记录。科学失败和阻断任务可保留在快照中供解释，但不能贡献接纳的证据。

快照包含分段内容摘要，用于发现损坏或普通编辑；摘要不是可信签名，也不证明科学结论真实。报告重建要求原注册表文件可用，且当前产物集合与快照一致：来源缺失、内容改变、同研究新增产物、审计撤销或收据丢失都会明确拒绝。快照本身不能修改或重新授予审计状态。

## 验证

新增测试覆盖只读预览与执行解析一致性、真实 h5ad 导入、严格恢复、不同工作目录调用、恢复不重导入外部文件、来源篡改、审计撤销、缺失计算收据、快照原子写失败、旧研究目录保护及清理阻断退出码。

2026-09-19 最终全量回归：**282 passed，0 failed，0 skipped**，耗时 198.07 秒；419 条警告涉及科学库弃用提示、数值边界及合成数据测试。命令为 `.venv\Scripts\python.exe -m pytest tests -q -p no:cacheprovider --basetemp=.p3_final_validation --junitxml=p3_final_validation.xml`，本机机器可读结果保存在 `p3_final_validation.xml`。`git diff --check` 通过。

安装验收使用离线构建的 wheel，在独立虚拟环境中安装，并从源码目录之外调用 CLI；科学依赖复用本机已有安装，因此这不代表全新机器的依赖解析验收。600 个细胞、300 个基因的合成数据运行包含知识检索，共完成 10 个任务，生成 7 个证据节点及 7 个结论；严格恢复复用全部 10 个任务且无重新计算，快照重建报告与原报告逐字节一致。

端到端验收还修复了知识检索的两个问题：发现模式现在与先验引导模式一样产生必需的 `epistemic_tagging_check`，并直接核对持久化报告的模式与标签，不能只靠任务返回的指标通过；无检索结果表的占位类别改为 `No_evidence`，避免 CSV 将字符串 `None` 解析为空值而误报统计审计失败。新增两种模式的真实产物与伪造指标回归测试。
