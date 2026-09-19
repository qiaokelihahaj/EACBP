# 架构检查与修复记录

## 2026-09-16：多产物任务事务

新增 TaskArtifactTransaction：每次执行/重试使用独立目录，输入从主注册表读取，输出仅在私有索引中登记。任务成功时统一校验 payload 哈希，在主索引一次原子替换中发布所有元数据与 TaskResult 恢复记录。发布前失败不占用正式 URI，发布后任务日志未保存也能恢复。保留原有审计与证据准入门禁。

测试覆盖真实子进程在首个输出后异常退出、发布后异常退出、索引写入失败时历史数据完整、首次任务日志保存失败后无需重算。未提交目录保留、不自动回收；外部工具副作用不在事务范围内，旧版本部分输出不能自动迁移。

本轮最终完整回归：**156 passed，131 warnings，95.80 秒**。差异格式检查与 wheel 构建通过；JUnit 结果在 `outputs/review_20260916/transaction-tests.xml`。

## 2026-09-16：并发、崩溃恢复与结论语义补充

修复同一 Registry 实例读取可能重置未提交注册信息的问题，增加可重入锁并返回元数据副本。增加 computed/audited 两阶段任务日志，只有最终日志成功后才准入状态、规划和证据；审计异常或日志提交失败可以重新审计已保存输出。

重复全量回归还复现 Windows 索引读句柄阻止另一进程原子替换的问题；索引初始化与刷新现与写入共用进程锁，防止该间歇性提交失败。

统计结论要求与支持的统计证据摘要一致；解释级的已验证知识还须与该证据生物学上下文相符，避免借用无关基因显著性。自由改写暂不自动准入。

新增测试覆盖实例内并发、四进程八产物提交、进程崩溃释放研究锁、审计/日志故障恢复及无关结论拒绝。逐项验收与未验证边界见 `outputs/review_20260916/修复验收清单.md`。

最终完整回归：**152 passed，131 warnings，84.05 秒**；追加存储/并发回归 11 passed。`pip check` 与 `git diff --check` 通过；隔离构建 wheel 成功。测试日志为 `outputs/review_20260916/acceptance-tests.xml`。警告类别与前轮一致，未屏蔽。

## 2026-09-16：继续修复的实际结果

本节取代下方历史记录中的“未实现”状态。本轮在子代理再次因额度停止后，由主代理完成算法实现和集成验证。

- **实际算法已接入**：`harmonypy_v1` 调用 Harmony；`scanpy_leiden_umap_v1` 调用真实邻接图、Leiden 与 UMAP；`scanpy_dpt_v1` 调用 diffusion map/DPT，并支持真实 PAGA。真实数据默认采用 standard 配置，明确的模拟输入采用 baseline；也可显式选择配置。
- **CellRank 已接入并验证**：`cellrank_fate_v1` 使用 PseudotimeKernel 或 PrecomputedKernel 与 GPCCA 估计器计算命运概率。支持 DPT→CellRank 的完整 DAG；要求用户明确指定终末细胞。已用已知解析答案的转移矩阵验证概率数值。没有实现从原始 spliced/unspliced 自动估计 RNA velocity。
- **STARsolo 执行分支已实现**：多 lane 的 cDNA/条形码顺序、chemistry 参数、参考和输入哈希、独立运行目录、Matrix Market 导入、非零退出处理。只读取本次输出；损坏或缺失 filtered 矩阵不再自动替换为 raw。命令行入口支持 STAR 参数。
- **规划不再无条件套用 Microglia 流程**：未指定细胞类型时分析全部细胞；找不到指定类型时失败，不替换为数量最多的类型。审计后根据实际条件和供体数裁剪统计及其下游分支；没有明确根细胞时省略 DPT 并记录原因。报告包含规划决策。
- **扰动完整性修复**：DEG 表不再生成 100 个假细胞；无明确对比时不构造疾病签名；药物参考与推断网络标明来源、验证状态；任意 M1/M3 聚类编号不再被当作疾病状态。遗传扰动前后使用同一个 PCA 基，修复把不同坐标基相加的错误。
- **恢复与输入溯源**：严格恢复检查源码与安装依赖指纹。补齐空间模拟工厂缺失的模拟标记。新增通用 `scripts/run_standard_study.py`，从真实存在的 h5ad 和显式研究设计生成报告，不依赖 Kat8 模板。

### 本轮验证

- 完整回归：**145 passed，80.56 秒**，包含真实 Harmony、Leiden/UMAP、DPT/PAGA、CellRank、审计与恢复执行。
- STARsolo 的 4 项测试使用受控的 subprocess 替身和 Matrix Market 文件验证命令、导入、失败退出和输入篡改；**没有运行真实测序对齐**。
- 扰动完整性及方法相关 10 项测试通过。完整回归仍有 131 条警告，主要是 AnnData 索引转换、Leiden 后端建议和近常量数据的 SciPy 精度提示。
- 环境：Windows / Python 3.13.13，AnnData 0.12.19、Scanpy 1.12.4、HarmonyPy 0.0.10、CellRank 2.3.2、Pandas 2.3.3、SciPy 1.16.3。完整安装快照见 `requirements-tested-windows-py313.txt`。
- CellRank 与 AnnData 0.13 的内部接口不兼容，因此 `fate` 依赖约束为 AnnData `<0.13`；HarmonyPy 2.0 需要本机 C++ 构建工具，因此使用已实测的 0.0.10 实现。依赖检查通过。

### 仍需外部验证的边界

真实 FASTQ 对齐需要实际测序输入、参考与可执行程序，本轮未执行。标准算法调用通过不代表其适合每一种生物学实验设计；DPT 稳定性未评估时明确不产出稳定轨迹结论。知识库仍是本地未验证内容；CellRank 概率依赖所选转移模型与终末细胞。没有把这些范围宣称为已经验证。

算法接口依据：[Scanpy Leiden](https://scanpy.readthedocs.io/en/latest/api/generated/scanpy.tl.leiden.html)、[Scanpy PAGA](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.paga.html)、[HarmonyPy](https://github.com/slowkow/harmonypy)、[CellRank GPCCA](https://cellrank.readthedocs.io/en/latest/api/_autosummary/estimators/cellrank.estimators.GPCCA.html)。

## 2026-09-15 历史记录（以下为当时状态）

本轮由三个 Luna Max 子代理分工修改输入、方法和存储，主代理负责执行治理、证据、报告与集成验证。子代理后期因额度限制停止，主代理接手完成回归修复。改动尚未提交 Git。

## 已修复

| 问题 | 当前行为 |
|---|---|
| 缺少 FASTQ 或工具时静默生成数据 | 默认 real 模式明确失败；只有显式 demo 生成合成数据，来源沿血缘传播至报告 |
| 输入遗漏 lanes、缺少设计信息、复用旧结果 | 收集所有配对 lanes，要求 donor/condition 等元数据，运行目录隔离；记录命令与 FASTQ/参考内容哈希，检查执行期间变化 |
| 审计失败仍继续生成结论 | 所有输出均接受计算审计；失败任务不贡献证据，依赖任务被阻断 |
| 非显著结果、空支持或模板疾病结论 | 按实际结果提取证据；统计推断要求已审计的显著性证据；拒绝空或不存在的支持 ID |
| 本地知识伪装显著实验发现 | 移除虚构 p 值，标为未验证背景假设；不自动提升为机制结论 |
| 文件篡改与重启丢失元数据 | 持久化索引/血缘，读取校验 SHA-256，原子发布、进程锁与 URI/循环血缘检查 |
| AnnData 数据槽丢失 | 保存 layers、raw、obsp、varm、varp 与嵌套元数据；修复 raw/varm 子集和新版 AnnData None layer 兼容 |
| 序列化掩盖错误 | 不静默降级损坏的 h5ad；无依赖时使用明确 NPZ 格式与禁止 pickle 的结构化编码 |
| 简化方法冒用第三方算法名称 | 使用真实方法 ID，保留显式迁移别名；K-means、批次均值中心化和根距离排序均不再宣称为 Leiden/Harmony/CellRank |
| counts 丢失与细胞伪重复 | 保留 counts；供体聚合后检验，细胞级回退明确标为探索性；拒绝不适用的配对设计 |
| 注释泄漏、空轨迹崩溃、空间二次方内存 | 不读取 ground truth 作为预测；空结果有稳定表结构；空间邻接使用稀疏近邻，扰动计算设规模边界 |
| 未知药物必然得到“逆转”结果 | 无已知签名时要求显式 drug_signature，不再合成反向疾病签名 |
| 无恢复和跨运行状态混用 | 按依赖排序、有限重试、任务日志、输入/输出哈希匹配后恢复并重新审计；每次调用清空内存研究状态 |

## 验证记录

- 本地 Windows / Python 3.13.13；安装真实 AnnData 0.13.3.post0、Scanpy 1.12.4 和 h5py 3.16.0 后运行。
- 完整回归：128 passed（48.15 秒）。随后新增 1 项未知药物回归并修改 FASTQ 内容溯源，相关输入/方法测试：12 passed。当前收集总数为 129；并非最后一次完整运行了 129 项。
- 完整测试存在 53 条警告，主要是 AnnData 字符串索引转换，以及一个近常量数据的 SciPy 精度警告；没有将警告屏蔽。
- `pip check`：无损坏依赖；wheel 构建通过；`git diff --check` 无空白错误。
- Kat8 脚本显式 demo 完整运行成功，报告在 `outputs/runs/Kat8_P12_cKO_SingleCell_Study_001/verified_bio_20260915/reports/Kat8_cKO_Study_Report.md`。
- CI 配置覆盖 Windows/Linux、Python 3.11/3.12/3.13；远端 CI 尚未执行。

## 尚需真实数据验证的边界

没有运行真实 FASTQ 对齐或 STARsolo；后者明确未实现。真实对齐的命令、参考兼容和生物学结果仍需真实输入验证。当前流程仍是固定能力组成的本地研究原型，不是通用自适应研究系统。

批次校正、聚类、轨迹、标记字典、空间分析和扰动模型仍为简化实现。内置药物/知识参考是本地人工数据，不能据此宣称疗效、实验因果关系或已验证机制；探索性统计也不替代适合实验设计的正式模型。哈希校验保护相对于本地索引的一致性，不提供外部可信签名。

原始检查报告与重现实验保存在 `outputs/review_20260915/`。历史规划文档已加说明，当前行为以 README 和可执行代码为准。
