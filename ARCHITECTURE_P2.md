# P2：能力扩展、运行维护、执行观测与多目标分析

本轮在 P1 的类型化配置、任务执行器和持久化审计边界上继续改进。保留本地同步执行方式和单目标任务兼容性。

## 能力声明

`CapabilityDescriptor` 汇集方法 ID、参数模型、输入/输出类型、必需审计项及可选的规划、独立审计和证据提取函数。注册表在计算签名前规范化任务参数，执行前检查输入和预期输出类型，执行后检查实际输出类型。必需审计项进入原有审计覆盖检查，失败结果不能生成证据。

六个高级扩展的声明集中在 `eacbp/capabilities/advanced_descriptors.py`，规划器从中读取方法、输出和审计要求。历史能力自动生成兼容声明；其旧类型元数据尚未全部迁移为强制类型约束，高级方法的细粒度科学参数约束仍由现有实现与审计器检查。

新扩展注册一个带 descriptor 的 `BaseCapability` 即可通过 `analysis_extensions` 请求，无须修改主编排循环。扩展函数约定：

- `plan_factory(manifest=..., parameters=..., tasks=...)` 返回该能力的一个或多个 `TaskContract`。
- `validator(contract, result, registry)` 返回对应任务的 `ValidationReport`，应独立读取持久化结果验证。
- `evidence_extractor(contract, result, report, registry)` 返回 `EvidenceNode` 列表，仅在审计通过后调用。证据必须引用当前任务的已注册输出，不能借用其他任务的来源或重复证据 ID。

参数模型可采用 Pydantic；应允许 `study_id`、`target_branch` 等编排上下文字段，或在模型中明确声明。完整的独立扩展例子见 `tests/test_capability_descriptors.py`。高级内置扩展的规划工厂提供输出骨架，仍由高级规划器负责其前后处理依赖。

## 执行事件

每次 `run_study()` 调用创建独立的 `run_id` 与 JSONL 日志，恢复运行也使用新日志，不覆盖上一次记录。返回摘要增加 `event_log` 和 `observability_warnings`。

事件覆盖规划、任务开始、方法选择、每次执行尝试、重试原因、产物提交、恢复命中、审计、依赖阻断和运行结束。任务与尝试耗时使用单调时钟，事件时间使用 UTC。日志逐条刷新；事件写入失败会进入返回摘要的警告列表，不改变科学审计结论。

```powershell
.\.venv\Scripts\python.exe scripts/inspect_run.py <event_log路径>
```

读取器允许忽略崩溃时最后一条不完整记录，但完整记录损坏、序号不连续或混入其他运行的事件会报错。事件日志用于诊断，任务检查点和审计收据仍是恢复与证据接纳的依据。

## 多目标观测隔离

共享预处理的观测和各目标细胞群的观测分别保存。某一分支的细胞数、批次等指标不会改变兄弟分支的状态；方法配置仍不能由返回指标覆盖或注入。目标和检验范围进入任务、结果和证据上下文。

`target_parameters` 可按目标细胞类型配置各能力参数，避免将同一个轨迹根细胞强加给不同细胞群。多目标报告保留各分支的分析范围，不将各自校正的显著性视作跨分支统一校正。

当 `manifest.biological_design.target_cell_types` 包含多个目标时，预处理与注释共享；子集、丰度、DEG、轨迹及其消费者使用独立任务 ID 和产物 URI。目标名形成的标识若冲突，会增加确定性哈希。单目标保持原有 URI。一个分支失败会阻断其下游，其他分支仍可完成；运行总状态仍按已有规则标记失败，避免把不完整研究报告为成功。

```python
manifest.biological_design.target_cell_types = ["Microglia", "Neurons"]
summary = orchestrator.run_study(manifest, {
    "method_profile": "standard",
    "target_parameters": {
        "Microglia": {"trajectory_inference": {"root_cell_id": "microglia_cell_1"}},
        "Neurons": {"trajectory_inference": {"root_cell_id": "neuron_cell_1"}},
    },
})
```

根细胞必须真实存在于相应子集中；缺少明确根细胞时不自动发明 DPT 根。启用多目标 CellRank 时，每个目标都需提供 `fate_mapping.terminal_states`。`fdr_family` 用于标识目标范围，实际多重检验校正仍分别发生在各项分析内，不提供跨分析或跨目标的联合 FDR 保证。

## 维护边界

清理仅面向注册表根目录中能够证明已废弃的未提交事务，不按 `temp*` 名称批量删除工作区目录。已提交产物的实际内容可能位于 `_transactions`，必须保留。

新事务使用操作系统锁标识活动状态；任务失败、完成或进程退出后释放。维护流程使用同一父注册表锁协调扫描与提交，并在执行删除前重新校验引用和路径。默认行为是预览；实际删除需要显式选择。

```powershell
# 默认只预览，保留最近 30 天内更新的事务
.\.venv\Scripts\python.exe scripts/cleanup_artifacts.py <artifact_root> --older-than-days 30
# 明确执行，重新扫描并校验当前状态
.\.venv\Scripts\python.exe scripts/cleanup_artifacts.py <artifact_root> --older-than-days 30 --apply
```

Python 接口 `plan_cleanup()` / `cleanup_artifacts(..., apply=True, plan=preview)` 支持绑定预览快照：索引、日志或事务内容变化会使旧计划失效，本次不删除。缺少或损坏的父索引、未知目录内容、符号链接、旧版无租约标记事务及仍有引用的事务均保守保留。输出包含候选数量、字节数和跳过理由。

## 验证

2026-09-19 全量验收：**256 passed，416 warnings，180.39 秒，无失败、无跳过**。相较 P1 的 236 项增加 20 项测试；额外 2 条警告来自新增多目标夹具的 AnnData 索引字符串转换，其余为既有科学库警告。

专项验证覆盖声明式扩展接入、参数及输出类型拒绝、审计失败不接纳证据、跨任务证据拒绝、双目标真实子集与 DEG、分支失败隔离与恢复、目标根与终末状态配置，以及事务保留期、活动租约、进程崩溃、提交引用、损坏索引和过期预览。

全量验收命令及结果文件：

```powershell
.\.venv\Scripts\python.exe -m pytest -q --basetemp=.pytest_temp_p2_final --junitxml=.pytest_cache/p2_final.xml
```

完整输出记录在 `.pytest_cache/p2_final.log`。本轮不清理现有研究数据；清理删除测试仅在测试夹具目录执行。未新增真实 FASTQ 比对、外部 CellBender 作业、远程 CI 或大规模性能基准。
