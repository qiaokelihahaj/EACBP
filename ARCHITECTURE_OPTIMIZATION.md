# 架构优化：运行配置、任务执行与审计状态

本轮范围是前三项高优先级改进：分离配置和观测状态、拆分执行循环、持久化产物审计状态。保留本地同步执行、现有能力实现、版本化产物和字典调用入口。多目标分支、分布式执行、统一能力描述、产物清理和证据快照不在本轮范围内。

## 运行时职责

| 对象 | 职责 |
|---|---|
| `RunConfig` | 调用者提供的运行选项与类型校验 |
| `PlanningContext` | 推导值、任务观测及规划决策 |
| `ExecutionState` | 运行中的状态与兼容旧 planner/router 的字典视图 |
| `TaskExecutor` | 执行单个任务、重试和方法回退，每次尝试使用私有事务 |
| `ResumeManager` | 签名生成、恢复结果检查、检查点内容构造 |
| `EvidenceAdmission` | 审计与暂存证据、规划变更，等待持久化成功后接纳 |
| `ScientificOrchestrator` | 按依赖调度并协调以上组件 |

配置和观测分开保存。任务返回的指标不能覆盖或注入方法选择、扩展启用等执行配置；原始指标仍可供诊断查看。批次和重复数等观测应由实际数据审计更新。调用者传入的嵌套配置与运行内部状态隔离。

字典入口和 `RunConfig` 入口保持相同的恢复语义。布尔控制参数要求真正的布尔值：例如 `resume="false"` 会在执行前报错，避免被 Python 的字符串真值判断误解为恢复请求。为兼容已有自定义能力，未知扩展字段仍允许传入；这不是完整的插件参数校验框架。

## 审计状态与持久化顺序

计算产物发布和科学审计是两个阶段。原始 `get()` 仍可读取计算产物；需要科学审计保证的调用方使用新增的受审计约束接口。

1. 在私有事务中计算，成功后原子发布产物及计算收据。
2. 保存计算检查点。
3. 开始审计，持久化 `pending` 状态，撤销同一任务上下文中的旧通过状态。
4. 执行独立审计并暂存证据及规划变更。
5. 保存审计阶段检查点，再持久化最终审计记录。
6. 成功后更新运行状态、计划和证据图。

审计记录关联任务签名、合同、产物哈希和审计器身份。未审计、待审计、被拒绝、上下文不匹配或文件完整性失败的产物不能通过受审计约束接口读取。旧注册表缺少审计记录时，不将已有计算产物默认视为通过审计。

恢复仍复用已提交的计算结果并重新审计。检查点或审计状态写入失败时，不应接纳新证据。严格复现模式仍包含源码指纹；因此升级代码后，旧运行可能需要新建运行目录，这不是跳过配置完整性检查的理由。

新增接口示例（`signature` 和 `resolved_contract` 必须对应已执行的任务）：

```python
from eacbp.schemas.runtime import RunConfig

summary = orchestrator.run_study(manifest, RunConfig(method_profile="standard"))

record = registry.get_audit_record(signature, contract=resolved_contract)
metadata, payload = registry.get_audited(
    output_uri,
    signature=signature,
    contract=resolved_contract,
    expected_auditor_fingerprint=expected_policy_fingerprint,
)
```

`get_audit_record()` 用于检查状态，本身不等于重新校验文件完整性；`get_audited()` 才是受审计约束的读取入口，并验证同任务的其他输出。调用者可以通过 `expected_auditor_version` 或 `expected_auditor_fingerprint` 要求指定审计策略。原始导入文件若没有对应计算与审计收据，仍使用普通读取接口。

## 边界

- 审计状态是本地持久化记录，不是外部签名或可信执行证明。
- 审计模块独立于计算模块，但仍在同一个 Python 进程运行。
- 自定义审计器的外部配置、运行时替换等若不能从源码识别，应提供明确的审计器身份。
- `_transactions` 中包含已提交产物，本轮不清理或移动既有数据。

## 验证

改动前关键回归：编排、任务事务、崩溃与并发、自适应规划、注册表注入，共 20 项通过。

2026-09-18 最终全量验收：**236 passed，414 warnings，165.54 秒**。

```powershell
.\.venv\Scripts\python.exe -m pytest -q --basetemp=.pytest_temp_arch_acceptance --junitxml=.pytest_cache/architecture_acceptance.xml
```

本地完整日志：`.pytest_cache/architecture_acceptance.log`；机器可读结果：`.pytest_cache/architecture_acceptance.xml`。

新增测试覆盖配置覆盖/注入、字典与类型化恢复入口、非布尔恢复参数拒绝、审计持久化失败后的状态隔离、旧索引兼容、审计上下文和审计器身份、同任务其他输出完整性、矛盾审计报告及重审撤销旧状态。原有事务、崩溃恢复及实际库分析测试一并通过。

本轮未另行执行真实 FASTQ 比对、外部 CellBender 作业、远程 CI 或大规模性能基准；本地套件通过不替代这些验证。
