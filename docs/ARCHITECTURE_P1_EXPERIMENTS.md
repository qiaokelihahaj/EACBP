# P1：逐拟合诊断、行级证据与 QC 敏感性

本轮在共享原生流程中落地三项改进，并在 F 盘公开数据上验证。代码留在现有仓库；输入派生文件、缓存、测试临时目录和实验结果均位于 `F:/EACBP-architecture-p1`。旧经典复现记录保留，本轮结果不替换其科学验收状态。

## 已落地的行为

1. **逐拟合凭据**：原生 advanced PyDESeq2 计划为主模型和供体留一任务增加注册 JSON 输出，与科学表一起进入计算提交、哈希和独立审计。保存输入计数、样本、真实设计矩阵、对比、有效参数、版本、实际离散度趋势、捕获的 Python 警告及失败折次。实际趋势取拟合状态 `dds.uns.disp_function_type`，不能用请求的 `fit_type` 替代。旧任务明确诊断不可用；完成或无警告均不证明收敛。子进程仅写入 stderr 的警告不在当前捕获范围。
2. **行级证据和反证**：native DEG / functional activity 证据绑定源制品 URI、唯一键和行摘要，逐项核对效应、SE、CI、FDR、估计状态、对比及群体。排序不改变定位，重复键、改数值或只改原生结论文字都被拒绝。抽取、快照写入、带 registry 的读取、报告渲染均校验源表；离线读取只校验结构。支持与反证显示同样的溯源和数值。汇总及旧证据明确 `unavailable`，不伪造行验证。
3. **QC 敏感性**：分析前冻结原始细胞全集、计数指标、原始文件 SHA、完整配置和主 3MAD / 较严 2MAD / 较宽 4MAD 的阈值。分支仅改变 `qc.min_genes`、`qc.max_mito_pct`，通过 `run_study` 分别完整重跑。比较核查实际过滤掩码、每个细胞的去留、共享细胞重聚类 ARI，以及完整基因集合的方向、FDR、缺失、不可估计状态和设计信息。主 snapshot / summary 不覆盖；派生比较单独保存。失败、预算耗尽、重复阈值、少于三个细胞和未准入结果不能冒充完成。

这些改进提高结果的可检查性，尚不构成独立生物学重复、因果证明或任意自然语言问题的完整理解。

## PBMC3k 全量实验

输入为 2700 × 32738 的公开原始矩阵，SHA-256 为 `89a96f1beaa2dd83a687666d3f19a4513ac27a2a2d12581fcd77afed7ea653a1`。三个方案都通过原生计算、审计和问题范围检查。参考标注、预期基因及固定评分门槛只在三个原生运行完成并校验后由独立评估进程读取。

| 方案 | 保留原始细胞 | 与主方案共享细胞 ARI | 参考标注巨核细胞保留 | 固定科学评分 |
| --- | ---: | ---: | ---: | --- |
| 主 3MAD | 2562 | — | 5 / 15 | 失败：每类保留率 |
| 较严 2MAD | 2450 | 0.8743 | 4 / 15 | 失败：总体及每类保留率 |
| 较宽 4MAD | 2637 | 0.8842 | 11 / 15 | 失败：每类保留率 |

较宽方案保留了更多参考稀有群体细胞，仍低于事先固定的每类 80% 门槛。三个方案均未完整通过科学评分；不调整默认阈值或评分门槛来制造通过。共享细胞 ARI 是同输入重聚类敏感性，表中的参考标注也是既有表达分析结果，不能把它们当作独立实验真值。本机版本的主方案参考 ARI 为 0.5318，不能宣称与旧 Linux 容器的 0.5633 位级一致。

当前冻结 v3 再次完整运行三个方案，结果与 v1 的固定科学评分相同。[数值一致性记录](F:/EACBP-architecture-p1/verification/pbmc3k-v3-replay-consistency.json)保留这次复跑比较；下图来自同数值的 v1 评分。

证据：[原生敏感性报告](F:/EACBP-architecture-p1/runs/pbmc3k-p1-v3/qc_sensitivity/report.md)、[固定参考科学评分](F:/EACBP-architecture-p1/verification/pbmc-qc-external-science-v3.json)、[运行凭据](F:/EACBP-architecture-p1/verification/pbmc3k-v3-experiment-receipt.json)。

![PBMC QC 群体保留率](F:/EACBP-architecture-p1/verification/pbmc-qc-retention.png)

## Kang 诊断链演练的范围及保留失败

本轮技术演练选择源文件已经记录的 B 细胞群，保留全部 2651 个 B 细胞、8 个配对供体和 15706 个基因。选择发生在结果计算前，依据观测群体而非预期基因或模型结果；因此这不替代旧七群体 Kang 全量复现。

最初输入仅有 `replicate`，新流程拒绝将通用重复列自动提升为供体，返回问题范围 `partial`。随后在独立输入派生步骤中显式记录 `replicate -> donor_id` 声明：原实验将八名个体分别分成刺激和对照 aliquot，[原论文实验设计](https://pmc.ncbi.nlm.nih.gov/articles/PMC5784859/) 与 [Pertpy 官方教程的 replicate / label 样本映射](https://pertpy.scverse.org/en/latest/tutorials/notebooks/milo.html) 提供依据。计数、原始细胞及完整特征空间逐项验证未变；输入文件 SHA 和角色声明一起记录。身份标签及声明本身仍不证明抽样独立性。该数据适配是明确的人工作业，不声称整个实验无需人工介入。

接入实验发现新 JSON 输出尚未进入能力类型/参数声明，且统一计划装配的多群体 URI 未同步。失败运行保持原状，修复后在新目录重新运行。既有原生集成测试捕获了这些问题，提醒流程改进必须经过声明、规划、执行、提交、审计、快照、报告整个路径，不能仅以能力函数测试成功为完成。

第二轮全部 30 次拟合完成，但单群体 DEG / 留一任务未带明确 B 细胞身份，严格问题完成检查仍返回 `partial`。修复规划时的身份传播后，在 v3 新目录完整重跑；三种 QC 下原生计算、审计及问题范围均通过。旧 v2 `partial` 保留。逐列精确比较三种 QC 的 DEG、留一逐折及汇总共 9 张表、518298 行，全部科学字段与 v2 相同，说明此次身份修复未改变统计数值。[比较凭据](F:/EACBP-architecture-p1/verification/kang-v2-v3-science-equality.json)记录每个源 URI、文件及 SHA。

| 方案 | 保留 B 细胞 | 共享细胞 ARI | 可配对估计基因 | 与主方案方向反转 | FDR 判定变化 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 主 3MAD | 2651 | — | — | — | — |
| 较严 2MAD | 2638 | 0.8770 | 13328 | 68 | 22 |
| 较宽 4MAD | 2651 | 1.0000 | 13333 | 0 | 0 |

三个方案各执行主 DEG、留一基线和 8 个完整供体删除折次，共 **30 次真实拟合**，每次保持 15706 个基因。独立验收逐项核对原源、派生源与保留子集的计数、细胞、特征、供体和条件标签。7 个预定 IFN 基因在全部 QC 下均为正效应、FDR < 0.05，并在每个方案的 8 个留一折次中保持此方向和显著性；从原计数独立重算的 CP10K 方向在八个供体中全部为正。这里只支持指定 B 细胞响应在这些扰动下稳定，不能推广到全部基因或独立队列。

主 DEG 表完整保留 15706 个输入特征，其中 13333 有原始 p 值、6675 有 FDR、2373 为不可估计。完整特征保留不等于每个基因都完成有效检验。供体敏感性制品仍记录 `scientific_robustness_claim_supported=false`；不能将分析范围完成提升为普遍稳健性结论。Kang 输入也没有可用的线粒体特征指标，本轮实际检验基因数阈值扰动，不能据此声称已检验线粒体 QC。

逐拟合诊断明确显示 **26 次 parametric、4 次 mean 回退**，四次回退都捕获了趋势拟合未收敛的 Python 警告：三个 QC 方案的 `patient_1015` 留一折次，以及较严方案的 `patient_1488` 折次。完成不等于收敛；较宽方案使用相同细胞集合，也不增加独立验证。

证据：[主原生报告](F:/EACBP-architecture-p1/runs/kang-b-declared-donor-p1-v3/report.md)、[QC 比较报告](F:/EACBP-architecture-p1/runs/kang-b-declared-donor-p1-v3/qc_sensitivity/report.md)、[逐拟合与独立科学核对](F:/EACBP-architecture-p1/verification/kang-b-diagnostic-verification-v3.json)、[运行凭据](F:/EACBP-architecture-p1/verification/kang-b-declared-donor-v3-experiment-receipt.json)。

只读复核还发现真实原生表将目标放在嵌套 `inference_contract`，早期行校验仅检查 flat 参数。最后补上从实际源行、嵌套源契约和任务身份的核对，并拒绝把证据自己添加的未绑定指标提升为来源身份。该修复在计算结束后应用，没有改写 v3 原生文件；新校验器重新验证六个完成快照，并在真实 Kang 节点上拒绝 12 次内存篡改。[实际源验证](F:/EACBP-architecture-p1/verification/real-native-target-tamper-verification.json)与[85 项定向回归](F:/EACBP-architecture-p1/verification/row-target-scope-tests-v4.log)保留完整记录。最终源码校验环境与原生计算时的源码指纹分别记录，不冒充同一源码执行。

## 验证与范围

最终源码核心 CI **775 passed**、标准流程 **154 passed**、高级统计 **71 passed**，无意外跳过；两个最终 CI 执行期间源码和包环境均保持不变。核心和标准分别保留一个既有显式 deselect；定向测试与这些套件重叠，不重复累计测试数量。

主要记录：[最终核心日志](F:/EACBP-architecture-p1/verification/core-p1-v4.log)、[最终标准/高级统计验证](F:/EACBP-architecture-p1/verification/capability-suites-v4.json)、[此前 38 项整合集成检查](F:/EACBP-architecture-p1/verification/integration-fix-v2.log)。

本轮是本机 Windows 原生运行；输入 SHA 与源码/包环境前后检查只提供本机可复核性，不宣称旧冻结容器验收、跨版本完全一致或远程执行证明。主 study export 不包含各个 QC 分支的独立 registry，需要保留整个 `qc_sensitivity` 目录及主运行。

## 由实验推动的下一步

- 将低 RNA 被排除细胞保留为可检查的候选集合，报告其群体相关损失和标记特征，区分阈值异常、技术损伤证据和可能的稀有群体。任何默认 QC 改动必须先在另一公开数据集及冻结标准上验证；不能用这次参考评分直接选择默认分支。
- 增加来源绑定的通用输入角色声明及验证，处理 `replicate` 等歧义列；需要提供实验依据，不能恢复通用别名猜测。
- 将 QC 派生执行、断点恢复、独立 export 和显式 QC 问题范围统一成可验证的研究套件。当前主完成状态与派生敏感性状态分开，QC 敏感性未加入主要求矩阵。
- 统一诊断制品的规划阶段：目前仅手动 `method_override` 而未开启 advanced 的 PyDESeq2 任务仍无自动凭据；将来以拟合凭据为输入的扩展还需要同步分支 producer 映射。失败主拟合的独立 attempt 日志也需与科学准入区分。

完整实验索引由 [F 盘验收记录](F:/EACBP-architecture-p1/verification/p1-acceptance.json)给出，包含六个原生运行、科学评分、逐拟合凭据、最终 CI 和十四份旧记录的未变核查。旧基准失败不因本轮工程改善而变成通过。
