# 科学研究价值验证：Kang 公共数据

**本轮结论：Kang 受控分析、donor 敏感性结果和自动化复核验收通过。** 两个 cell type 的全模型、八次 donor 剔除模型、独立 audit 和已完成 reference fits 的最终逐基因重核均通过。敏感性判读规则是在已查看全模型后、LOO 运行前冻结的，并非完整研究的事前预注册。验收范围限于该公开数据与这套规则；新发现、外部生物学验证、人工时间节省和用户实际课题仍待验证。

本轮外部数据和验收产物位于 `F:\EACBP-public-validation\research-value-20261001`。最终镜像、focused 测试、运行记录与 QA JSON 都保留在该目录下。

汇总证据索引见 [research-value-acceptance.json](F:/EACBP-public-validation/research-value-20261001/runs/research-value-acceptance.json)。其中分别记录方法验收通过和仍待领域复核的科学结论。

## 验收目标与判读约定

本轮评估的是：在保持统计设计、输出完整性和结果可复核的前提下，EACBP 能否把单细胞数据处理成可审查的候选结果，并降低机器处理成本。当前没有可信的人工操作计时，因此 `human_active_minutes=null`、`human_savings_claim_supported=false`；不把机器运行时间换算成人时，也不声称节省人工时间。

两个 Kang 2018 子集各保留全部 15,706 个基因和 8 个配对 donor。全模型包含 16 个 donor-condition pseudobulk 样本；每次完整 donor 剔除后保留 7 对、14 个样本。分析使用既定设计 `~condition + donor` 和对比 `stim vs ctrl`。每个模型均保留完整基因表；合法的不可估计基因及 NA 值仍保留。

LOO 是对单个 donor 影响的**选择后敏感性分析**，不是独立验证。候选按全模型 `q < 0.05` 且 `abs(log2 fold change) >= 1` 筛选；方向稳定要求 8/8 个 LOO fold 均有有限、非零且与全模型同方向的效应。q 值跨阈值变化单独描述，不作为假阳性裁决。7 个预设 ISG marker（ISG15、IFI6、LY6E、MX1、IFIT1、IFIT2、ISG20）单列为已知正例基准；其余候选仍属探索性结果，尚未做新颖性审查。FDR 家族按各自 cell type 的全基因 `stim vs ctrl` 对比独立计算，不据此声称 cell-type 特异性。

## 已通过的自动化验收

| 项目 | 结果 | 证据 |
|---|---|---|
| 最终容器镜像 | `eacbp:research-value-cpu`；`sha256:c5ad5cadcacc342249d4c975b284199eacc7c5b4590637ebfeecb34fa2236a60`；395,414,014 bytes | `references\docker-image-audit-accepted.json`；Linux amd64、Python 3.12.14、UID 10001 |
| 镜像源码身份 | 镜像安装的 111 个 Python 文件与当前源码逐文件匹配 | `references\source-identity-audit-accepted.json`，状态 passed |
| Windows focused tests | 112 passed、0 skipped，44.74 秒 | `runs\research-audit-final-focused-rerun.log/.xml` |
| Linux focused tests | 112 passed、0 skipped，23.92 秒 | `runs\research-linux-audit-accepted-focused.log/.xml` |
| CD14+ Monocytes 最终分析 | overall passed；9 个模型（全模型 + 8 个 LOO），8/8 LOO 成功、0 失败；每模型 15,706 个基因；ScientificAuditor passed；356.559 秒 | `runs\cd14-audit-final\run_summary.json`、`audit.json`、结果表和卡片 |
| CD4 T cells 最终分析 | overall passed；9 个模型（全模型 + 8 个 LOO），8/8 LOO 成功、0 失败；每模型 15,706 个基因；ScientificAuditor passed；362.908 秒 | `runs\cd4-audit-final\run_summary.json`、`audit.json`、结果表和卡片 |
| 最终 reference 重核 | passed；两个 cell type 各 9 个既有模型结果均通过比较。每个模型的 15,706 个基因、5 项统计量及 NA/Inf mask 一致；输入 SHA、锁定环境、样本数、donor 数、design rank 与 residual degrees of freedom 也核对通过 | `runs\final-reference-recheck.json`，SHA256 `8ea190c9c85fa25910beb013ead6afd00b1292d3a36e5f8a3f3932c5bbb4f925` |

最终 reference 重核复用了同一输入 SHA 和锁定环境下已经完成的 18 个独立 PyDESeq2 reference fits，没有新拟合 reference。该核对确认数值实现一致，不等同于独立生物学验证。

## 敏感性结果与具体例子

最终结果卡排除 7 个已知 marker 后，CD14+ Monocytes 有 2,153 个探索性候选、CD4 T cells 有 470 个。分别有 381 和 79 个候选在 LOO 中出现 q 阈值状态变化；两组均无方向反转，也无效应覆盖不完整的候选。所有 7 个已知 marker 在两个 cell type 中都呈预设正方向，并且每个 marker 都有 8/8 有限效应、8/8 同方向效应和 8/8 个 LOO q 值低于 0.05。这验证了预设正例基准，不是新颖性或发现检验。

全基因表的有限效应覆盖另行记录：CD14 中 12,868 个基因具备 8/8 有限 LOO 效应，2,838 个不完整；CD4 中对应为 14,029 和 1,677 个。合法 NA 原样保留。拟合完成的检查核对每个模型的完整基因行表与设计，而非要求每个基因均可估计。

两张探索性结果卡展示了 LOO 对候选排序和解释的实际用途：效应方向可以稳定，但全模型的 q 阈值通过情况仍可能受个别 donor 影响。

样本数减少、各 fold 的独立过滤和检验功效变化也会改变 q 值。阈值变化不能单独归因于 donor 异常；下表的影响指标仅表示剔除该 donor 后的最大绝对 LFC 变化。

| Cell type / gene | 全模型 log2FC（95% CI）；q | LOO 效应方向 | LOO q < 0.05 | 最大绝对 LFC 变化；对应的留出供体 | 解释 |
|---|---|---:|---:|---|---|
| CD14+ Monocytes / **ENGASE** | -2.300（-4.121, -0.479）；q=0.03895 | 有限 8/8，同方向 8/8 | 2/8；6 位 donor 的 q 状态变化 | 0.80745；`patient_1015` | 方向在全部 donor 剔除中一致，但显著性阈值状态较敏感；q 变化不是假阳性裁决 |
| CD4 T cells / **PYCR1** | 2.356（0.964, 3.748）；q=0.007797 | 有限 8/8，同方向 8/8 | 7/8；1 位 donor 的 q 状态变化 | 0.50647；`patient_1488`（q 状态变化来自 `patient_1244`） | 效应方向一致，个别 donor 改变 q 阈值判断；仍需领域解释与外部验证 |

ENGASE 和 PYCR1 均标为 `post_selection_exploratory_sensitivity`，结果卡中 `novelty_claimed=false`、`independent_validation=false`。这两个例子来自分别分析的 cell type；不能用两个列表差异宣称 cell-type 特异性，也不能据此声称发现了新基因。

## 机器时间与适用范围

一次合成局部性能对照覆盖 3,000 个基因和 8 个 LOO fold。本次测的是 **LOO 摘要生成与摘要一致性核验两个 helper**：此前 19.699606 秒，本次 3.765992 秒，约 5.23 倍；结果 dtype 和数值逐位一致，audit errors 为空。分项本次为 1.1788855 秒与 2.5871064 秒，记录见 `references\loo-performance-after-final.json`。这不是完整 auditor、完整 workflow 或真实 Kang 全流程的速度对照，不能外推成人工节省。

最终镜像的 cached-wheel rebuild 记录为 105.1512544 秒，见 `references\docker-build-audit-accepted-timing.json`。冷科学依赖安装没有计时；镜像构建时间与 Kang 运行时间分开记录，不合并为端到端提速。

此前的 `runs\cd14-final`、`runs\cd4-final` 是 audit 失败的诊断运行，保留作历史记录但不作为最终结果。首轮停止运行的计数不一致记录在 `references\stopped-first-runs.json`：一次 donor 剔除后，旧输出错误沿用了全数据 8 对、16 样本计数，实际应为 7 对、14 样本。当前结果已通过更严格的设计与审计核对；保留历史失败是为了呈现修正经过。

## 尚未验证的科学价值

这一结果支持在 Kang 公共数据上复现固定设计、检查 donor 敏感性，并生成可追溯的全基因候选卡。它没有回答候选是否新颖、机制是否合理、能否迁移到其他队列，或 EACBP 是否减少人工分析时间。

下一步应在用户的实际课题中先确认科学问题、独立单位、estimand、对比和多重检验家族，再决定现有分析组合是否适用；对课题候选做领域复核，适用时进行独立队列验证；并为同一科学问题设计人工基线计时。ENGASE、PYCR1 是本轮方法示例，实际课题的候选应由课题问题与数据决定。完成这些工作之前，结论保持在“Kang 方法与敏感性自动化验收通过”，不升级为新发现或人工效率提升。
