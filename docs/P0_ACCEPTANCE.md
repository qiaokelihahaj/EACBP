# 下一阶段修正与验收

本轮依据《规划 EACBP 开发》先落实三个 P0，保留现有编排、事务、恢复、独立审计与证据准入机制。工作区原有修改不在本轮重新归因。

## P0 的验收边界

| 主题 | 必须验证的行为 |
| --- | --- |
| 安装与 CI | core 与可选科学依赖分开；各能力作业显式安装依赖；指定测试未收集、跳过或失败都不能算通过 |
| Contrast / alpha | 反向对比的效应、置信区间、标签一致；数值 contrast 必须明确效应定义；计算、审计、证据准入和报告使用相同 alpha |
| 稀疏 pseudobulk | 执行与审计均不展开细胞级完整 counts；校验非法 counts；聚合结果与稠密参考一致；资源预检在分配前生效 |

本地测试使用已有 Windows / Python 3.13 环境。它不能证明干净 Linux / Python 3.12 安装、远端 GitHub Actions、容器或 Slurm 已通过。发布锁文件应来自干净目标环境的成功解析与能力测试，不从 Windows 快照推导。

## P0-3 稀疏 pseudobulk RSS 探针

可在独立进程中重跑合成稀疏计数基准：

```powershell
.\.venv\Scripts\python.exe scripts/benchmark_pseudobulk_memory.py
```

2026-09-27，Windows / Python 3.13 默认参数为 200,000 个细胞、30,000 个基因、200,000 个非零计数。若将细胞级 int64 counts 完整展开，矩阵需要 48,000,000,000 bytes（48 GB，约 44.7 GiB）；实际 pseudobulk 表为 6×30,000。计数输出估算为 1,440,000 bytes，估算的聚合工作区为 13,920,000 bytes。2 ms RSS 采样从 168,243,200 bytes 基线观测到 206,475,264 bytes 峰值，增量 38,232,064 bytes。观测增量包含 pandas 分组及 Python 临时对象；估算值只表示计数聚合工作区，两者统计范围不同。

同一尺寸的独立审计输入重建可用 `.\.venv\Scripts\python.exe scripts/benchmark_pseudobulk_memory.py --operation audit` 验证；本次 RSS 基线为 168,615,936 bytes，采样峰值为 208,568,320 bytes，增量为 39,952,384 bytes。该路径校验原始 counts 并重建分组与设计矩阵，仍保持细胞级 counts 稀疏；它不是完整结果审计或模型拟合的内存测量。

此探针记录的是轮询得到的 WorkingSetSize，短于采样间隔的峰值可能漏记，因此不代表操作系统精确的进程高水位。稀疏回归测试命令 `.\.venv\Scripts\python.exe -m pytest tests/test_advanced_statistics_sparse.py --basetemp=.pytest-basetemp-sparse-20260927-03 -q` 完成，11 项通过。

## 本地验证记录

2026-09-27，在已有 Windows / Python 3.13 科学环境中：

- 全套回归：498 passed、1 skipped（缺少真实 CellBender smoke 报告），无失败；此前并发修改期间出现的通信恢复指纹失败，在源码稳定后通过。
- 随后补强非法 alpha 与缺失置信区间列的拒绝逻辑，最终运行 `scripts/verify_ci_capability_tests.py advanced-statistics`：44 passed、零跳过，包含稀疏、contrast、审计与端到端回归。
- wheel 构建、JavaScript 语法检查与修改文件 whitespace 检查通过。
- 仅含 core 依赖的干净虚拟环境已创建，但依赖下载约 17.7 kB/s，安装中止，未执行干净环境测试。远端 CI 和 Linux/Python 3.12 干净安装仍待运行；未生成发布依赖锁文件。

## 后续研究验收

下面是后续里程碑，不能用合成测试的通过替代：

1. 固定一套真实多供体 h5ad 与元数据，记录来源、许可、校验和、物种、供体、条件、批次、配对和原始 counts 层；当前没有为本轮提供这样的固定数据。
2. 在固定科学依赖环境中，比较 EACBP 与直接调用 PyDESeq2 的效应、区间和校正结果；记录数值容差及参考命令。
3. 验收全部不显著时正常完成、一个供体时拒绝供体级推断、完全混杂时拒绝不可识别模型、中断恢复不重复发布证据。
4. 对真实数据记录执行和审计峰值内存、耗时、细胞数、基因数、非零元素数、供体分组数及资源限制。

P1 的最小推断契约与结果状态、研究验收包、容器/Slurm，以及 P2 结果审查界面按上述验收结果继续推进。不显著结果只能说明当前数据未提供足够差异证据；实际等效需要预定界限及相应检验。
