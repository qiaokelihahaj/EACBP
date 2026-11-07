# P3：聚类执行收据与预先固定的稳定性实验

本轮沿用 PBMC3k 和 Baron human pancreas 的公开计数输入，检验原生聚类在固定扰动下的变化。沿用的原始输入分别位于 `F:/EACBP-autonomous-practice/datasets` 和 `F:/EACBP-architecture-p2/datasets`；本轮派生矩阵、运行目录、缓存、日志和验收记录位于 `F:/EACBP-architecture-p3`。本文件与源码位于工作仓库。P1、P2 的历史结果保留。

## 这次架构改进解决什么

旧聚类参数中的 `n_components` 是请求值。上游已有 `X_pca` 时，聚类直接使用它，邻接图会使用其中的全部维度，因此请求值不一定改变实际计算。进一步核对发现，现有 no-correction 稀疏路径用 `scipy.sparse.linalg.svds`，执行未中心化截断 SVD；密集路径和稀疏失败后的密集回退则中心化后调用 NumPy SVD。以前将这些路径统称 PCA，报告不足以说明实际行为。

本轮保留这些计算路径、随机种子和数值行为，补充上游表示与聚类执行收据。收据记录源制品哈希、细胞与特征顺序、计数来源、HVG 掩码、算法、请求与实际维度、邻居数、距离度量、图哈希、Leiden 分辨率、种子、实现与迭代次数、分区哈希，以及实际观测到的 donor/batch 列。未中心化 SVD 不再被报告成已执行中心化 PCA。

本轮 14 条原生收据均显示：clustering 请求 `n_components=20`，实际继承并使用上游表示全部 **30** 个维度；上游 no-correction 本身请求且实际输出 30 维。报告现在区分请求值与执行值，没有将这个现有行为误写成只用了 20 维。

`LeidenParameters`、`NoCorrectionParameters` 和 `HarmonyParameters` 限定各方法支持的参数。错误拼写和不适用参数会在执行前被拒绝。规划器仅移除自己为旧 k-means 注入的默认 `k_clusters`；用户显式给 real Leiden 传入这个参数仍会被拒绝。

`leiden_execution_integrity` 是 real Leiden 的强制审计。审计从输入表示重新构建邻接图、重新执行 Leiden，比较图与分区，核对输入计数、轴、参数、来源及汇总。簇名称整体重命名不影响分区等价性。继承的上游表示只核对算法路线、特征掩码、维度、哈希和父源，**没有再次执行上游 SVD**。审计复用同一第三方 Scanpy 库，检查执行重建；它不证明细胞类型真值、供体抽样独立性或优化已收敛。UMAP 为展示输出，不在重建范围内。

原生 Markdown 报告先要求唯一匹配的持久审计接纳，再读取制品、检查收据与合同及当前元数据的一致性。报告本身不重新拟合图。旧制品没有新收据时明确显示信息不可用。缺失 silhouette 也显示未计算，不再用数值 0 代替缺失。

## 固定的实验与实现边界

这是独立应用 API，未把额外聚类分支自动加入每次 `analyze`。调用者提供 count-only 自主探索配置与 manifest：

```python
from eacbp.application.clustering_sensitivity_execution import run_clustering_sensitivity

result = run_clustering_sensitivity(
    manifest=manifest,
    config=config,
    data=source_h5ad,
    output_dir=experiment_directory,
)
```

需要跨数据集统一冻结时，先调用 `prepare_clustering_sensitivity`，保存返回的文件哈希，再调用 `execute_clustering_sensitivity`。两者将计划、实际源文件、派生完整计数矩阵、实际配置、manifest 和运行环境绑定在第一次原生运行之前。执行直接调用 `run_study`，保留每个分支自己的 registry、持久审计、snapshot、summary 和 report。

v1 限于 standard 原生探索。主分析固定 resolution=1、n_neighbors=15、聚类 seed=42，上游 integration/PCA seed=42。扰动如下：

| 扰动 | 改变 | 保持 |
| :--- | :--- | :--- |
| resolution 0.5、2.0 | Leiden 分辨率 | 原计数、QC、归一化、表示算法及种子 |
| clustering seed 17 | 邻接图、Leiden、UMAP 的种子 | 上游表示种子仍为 42 |
| 原始细胞 80% 抽样 | QC 前从完整原始细胞抽样，seed=20261003，保留原顺序与全部基因 | 主分析 QC 阈值与科学配置；子集重新拟合 HVG 与表示 |
| 每个供体留出一次 | QC 前删除一个完整观测供体 | 固定 QC 与科学配置；子集重新拟合 |

供体留出仅在完整观测 donor 元数据含 2–4 个供体时计划全部留出。PBMC 缺失 donor，明确不可用；Baron 有四个供体，执行四次留出。样本、文库和技术批次标签不会自动变成供体。供体标签本身不证明抽样独立性。

PBMC 主 QC 固定 min_genes=400、max_mito_pct=5；Baron 固定 685、100。子集不重新用其分布选择新阈值。排除候选保持独立，不能回流到聚类。每个分支实际计数、完整特征轴、选择细胞及顺序、QC 保留细胞必须与运行前计划一致，全部原生任务与制品通过持久接纳后才可进入比较。失败、缺失或不可估计分支保留原状，没有补救阈值、获胜分支或默认结果替换。

## 结果应如何解释

比较使用主分析和分支共同保留的细胞，返回 ARI、采用 arithmetic 平均归一化的 NMI，以及最大化共同细胞 Jaccard 总和的 Hungarian 一对一簇匹配。每个主簇都保留完整细胞数、共同细胞数和覆盖率；匹配保留交集、并集和分支共同细胞数。没有匹配时得分为 0，单簇没有共同细胞时相关估计为 null；总体共同细胞少于 2 时 ARI/NMI 不可估计，返回 null。共同支持少于 20 个细胞会标记；单簇退化情况单独显示。

这些数值是依赖同一数据的扰动描述，不是独立生物验证。高总体 ARI 可能同时存在一个覆盖低或匹配差的小群体。LOO 中有意移除供体造成的覆盖下降必须与 QC 排除区分。此轮不规定“稳定性通过”阈值，也不根据外部标签选分辨率或改变默认流程。

外层 `comparison.status=completed` 只表示预定且可估计的分支已执行并完成比较。原生科学问题完成状态、外部重现科学准则、优化收敛和生物稳定性是不同判断。旁路比较不会补足原生问题合同中尚不支持的 robustness 能力。

团队此前已经见过 P2 外部评分，因此本轮只能声称冻结及原生执行未使用参考标签。随后同数据集参考比较属于探索性跟进，不能称为独立确认性盲验证。外部比较同时保留分支原始细胞/完整原始细胞覆盖、QC 保留/分支原始细胞保留和参考标签覆盖，避免混淆分母。

## 协议与证据入口

- [运行前冻结协议](F:/EACBP-architecture-p3/verification/p3-clustering-stability-protocol-v1.json)
- [运行前细胞名单](F:/EACBP-architecture-p3/verification/p3-branch-cell-masks-v1.json)
- [计数输入证明](F:/EACBP-architecture-p3/verification/p3-count-only-input-proof-v1.json)
- [原生运行前实际文件、配置与环境补充冻结](F:/EACBP-architecture-p3/verification/p3-native-freeze-supplement-v1.json)
- [独立重读 14 条原生运行并重算 12 组比较](F:/EACBP-architecture-p3/verification/p3-independent-native-stability-v4.json)
- [P3 对 P2 主分析科学结果零容差回归](F:/EACBP-architecture-p3/verification/p3-primary-science-regression-v1.json)
- [历史 20 份指定记录的字节保护检查](F:/EACBP-architecture-p3/verification/p3-historical-protection-final-v1.json)
- [外层 launcher 失败后的独立汇总恢复](F:/EACBP-architecture-p3/verification/p3-launcher-postprocessing-recovery-v1.json)
- [全部 14 条探索性外部科学评分](F:/EACBP-architecture-p3/verification/p3-exploratory-external-science-v1.json)及[紧凑逐分支摘要](F:/EACBP-architecture-p3/verification/p3-exploratory-external-science-v1.summary.log)
- [本轮完整验收索引](F:/EACBP-architecture-p3/verification/p3-acceptance.json)及[索引完整性检查](F:/EACBP-architecture-p3/verification/p3-index-integrity-final.json)
- [图形导出收据](F:/EACBP-architecture-p3/figures/figure-receipt-v1.json)及[视觉检查](F:/EACBP-architecture-p3/verification/p3-figure-visual-qa-v1.json)

协议固定 PBMC 5 条、Baron 9 条原生分析。实际文件、配置、manifest、命令与环境另有原生运行前补充收据；协议原文件保持不变。

官方方法说明：[Scanpy Leiden](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.leiden.html)、[SciPy svds](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.svds.html)、[scikit-learn ARI](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.adjusted_rand_score.html)。实际依赖版本以冻结环境为准，文档站点版本不替代本机运行收据。

## 本机实测结果

最终三条 CI lane 共 **1206 passed、0 failed、0 error、0 skipped**：core 920、standard 215、advanced-statistics 71。core 与 standard 各有一项显式 deselect，与 P2 的测试范围相同；“0 skipped”不表示所有可选测试都在本机运行。新增检查覆盖参数拒绝、执行收据、可重命名分区、持久接纳、有效重哈希后的篡改、缺失证据及完整分支流程。

PBMC 的 5 条和 Baron 的 9 条原生分析全部成功，原探索问题状态均为 complete，预定 12 组主分析对分支比较均完成。两套主分析对 P2 的独立回归核验通过：QC、归一化、integration、clustering 和 marker 的科学数组、细胞/特征顺序、计数 dtype 及表格均按零容差比较，没有数值变化。新执行收据没有以更换默认科学结果来制造改进。

独立校验 v4 已通过，逐一重读 14 条运行的持久接纳与全部制品，核对完整计数、细胞/特征顺序、QC 分区、供体角色、实际执行收据，并用独立实现重算全部 12 组全局指标、簇匹配与供体覆盖；369 份指定输入/原生/冻结文件在前后核验中保持相同哈希。它没有导入新聚类比较模块，指标采用 scikit-learn 与 SciPy；仍属于同数据、同第三方方法生态中的工程和描述性核验。参考标签与外部评分在该 gate 通过前没有打开。

历史保护检查另覆盖预先列出的 20 份 P1/P2/源记录，并全部保持字节一致；该范围不代表穷举了所有历史目录。独立验证器的 24 个合成病例和外部评分器的 8 个合成病例属于证据工具检查，未计入上述 1206 条仓库 CI 测试。

| 数据集与扰动 | 原始分支细胞 | QC 保留 / 共同细胞 | 共同细胞 ARI | 共同细胞 NMI |
| :--- | ---: | ---: | ---: | ---: |
| PBMC resolution 0.5 | 2700 | 2562 | 0.8377 | 0.8785 |
| PBMC resolution 2.0 | 2700 | 2562 | 0.7039 | 0.8419 |
| PBMC seed 17 | 2700 | 2562 | 0.9271 | 0.9267 |
| PBMC 原始细胞 80% | 2160 | 2046 | 0.8441 | 0.8787 |
| Baron resolution 0.5 | 8569 | 8569 | 0.8312 | 0.9126 |
| Baron resolution 2.0 | 8569 | 8569 | 0.8082 | 0.9348 |
| Baron seed 17 | 8569 | 8569 | 0.9343 | 0.9630 |
| Baron 原始细胞 80% | 6855 | 6855 | 0.9428 | 0.9551 |
| Baron 留出 GSM2230757 | 6632 | 6632 | 0.9428 | 0.9477 |
| Baron 留出 GSM2230758 | 6845 | 6845 | 0.9352 | 0.9497 |
| Baron 留出 GSM2230759 | 4964 | 4964 | 0.8906 | 0.9378 |
| Baron 留出 GSM2230760 | 7266 | 7266 | 0.9235 | 0.9415 |

以上比较分母是共同保留细胞。主分析 PBMC 保留 2562/2700；Baron 保留 8569/8569。原始 80% 抽样和供体留出造成的源细胞减少不能记为 QC 损失。改变分辨率时，ARI 下降部分体现簇拆分/合并；这些分数不能直接推导某个分辨率更正确。

![预先固定的聚类扰动结果](F:/EACBP-architecture-p3/figures/clustering-stability-v1.png)

图中间列显示共同细胞上的一对一匹配 Jaccard，右列显示共同细胞占完整主簇的覆盖率。星号标记共同支持少于 20 个细胞。高总体一致性与局部不稳定可以同时出现，应结合各簇覆盖、未匹配和完整分母审阅。

例如 Baron 留出 GSM2230758 后总体 ARI=0.9352，但主簇 `8` 的 326 个细胞仅剩 1 个，共同覆盖 0.31%，其一对一匹配交集为 0。这个分支没有 QC 排除，减少来自预定供体留出；总体分数不能据此宣称该群体跨供体稳定。匹配行存在也不表示相同群体：Hungarian 可分配零交集配对，实际交集/Jaccard 必须同时审阅。簇编号是当前运行的分区标识，没有被赋予细胞类型真值。

## 探索性外部科学评分

工程核验 14/14 通过；按历史冻结的完整科学准则，**1/14 通过、13/14 未全部通过**。两套 primary 均未通过完整科学准则。预先固定的 Baron resolution=0.5 分支通过旧门槛，外部 ARI=0.507114、NMI=0.774512，整体/各作者类别覆盖与七个 marker 描述也通过。它是一个被保留的探索性发现，没有被选作 winner，没有替换主分析或修改默认参数。此次供体留出分支均以原 resolution=1 为基准，**没有验证 resolution=0.5 的跨供体稳定性**。

| 数据集与分支 | 外部 ARI | 外部 NMI | 全部固定科学准则 | 未通过项 |
| :--- | ---: | ---: | :---: | :--- |
| PBMC primary | 0.5318 | 0.7260 | 未通过 | 每个参考类覆盖 |
| PBMC resolution 0.5 | 0.6159 | 0.7696 | 未通过 | 每个参考类覆盖 |
| PBMC resolution 2.0 | 0.3720 | 0.6643 | 未通过 | 每个参考类覆盖、ARI |
| PBMC seed 17 | 0.5198 | 0.7198 | 未通过 | 每个参考类覆盖 |
| PBMC 原始细胞 80% | 0.5783 | 0.7483 | 未通过 | 整体及每个参考类覆盖 |
| Baron primary | 0.4091 | 0.7155 | 未通过 | 参考聚类一致性 |
| Baron resolution 0.5 | 0.5071 | 0.7745 | **通过** | 无 |
| Baron resolution 2.0 | 0.3039 | 0.6808 | 未通过 | 参考聚类一致性 |
| Baron seed 17 | 0.3951 | 0.7095 | 未通过 | 参考聚类一致性 |
| Baron 原始细胞 80% | 0.4224 | 0.7249 | 未通过 | 整体/各作者类别覆盖、参考聚类一致性 |
| Baron 留出 GSM2230757 | 0.4959 | 0.7527 | 未通过 | 整体/各作者类别覆盖、参考聚类一致性 |
| Baron 留出 GSM2230758 | 0.4758 | 0.7384 | 未通过 | 整体/各作者类别覆盖、参考聚类一致性 |
| Baron 留出 GSM2230759 | 0.4261 | 0.7289 | 未通过 | 整体/各作者类别覆盖、参考聚类一致性 |
| Baron 留出 GSM2230760 | 0.4523 | 0.7210 | 未通过 | 各作者类别覆盖、参考聚类一致性 |

PBMC 的 reference label 仅覆盖原始 2700 个细胞中的 2638 个。未标注的 62 个细胞仍未标注，不补标签。原 source 全量分支的 Megakaryocytes 保留 **5/15**：10 个 low_gene_count 候选，不是供体遗漏。80% 分支的该类保留 **3/15**，其中源名单有意遗漏 5 个、真实 QC 候选 7 个；分支内 QC 保留率为 3/10，完整原参考覆盖为 3/15。两个分母都保留，排除候选没有回流。

Baron 九条分支都没有 QC 候选，subset 的类别覆盖下降全部来自预定源名单。80% 抽样按 floor 选 6855 个，其完整源覆盖 **6855/8569=0.7999766601**，略低于旧准则的 0.8，因此原整体覆盖检查仍失败；另有各类别覆盖及聚类一致性失败。这里不能将预定遗漏解释成 QC 质量损失，也不为通过而改旧准则。比如留出 GSM2230759 后作者 acinar 类从 958 个变为 115 个，完整源覆盖 12.00%，全部 115 个分支内细胞都通过 QC。

Baron 七个 published-marker 描述在九条分支均通过，但它们按作者标签计算 mean CP10K，而非独立细胞类型注释或推断检验。PBMC processed 标签及 Baron 作者标签同样依赖表达数据。一个分支达到这些固定门槛，不等于独立重现了原论文全部算法、图形或生物结论；也没有诊断候选细胞的损伤/活性或建立因果关系。

固定门槛仍沿用旧协议：PBMC 整体参考覆盖≥0.95、每类≥0.8、外部 ARI≥0.4、NMI≥0.65，以及至少六个受支持 marker 类别的原排名/纯度要求；Baron 整体/每作者类别覆盖≥0.8、ARI≥0.5、NMI≥0.6、七个 marker 配对全部满足原描述要求。两个数据集采用各自提前固定的准则，没有为了分支结果修改门槛。

完整验收索引绑定本轮源码、CI 收据、原生目录、冻结计划、独立校验、外部结果、失败版本与图形导出的文件哈希。审阅或备份需要一起保留原生 run 的 registry、持久审计、snapshot、制品、summary 和 report；本轮完整性保证限于索引与收据明确列出的本机文件，不代表外部可信证明或跨主机验证。

## 失败记录与资源边界

冻结的外层 launcher 在两套原生比较都完成之后，将字符串传给需要 `Path` 的哈希函数，构造最终汇总收据时退出码为 1。两份原始错误日志保留；冻结 launcher 没有重写，原生分析没有因此重新运行。独立补充收据仅重建源/环境一致性与原生完成字段，不能将 launcher 的退出码改成 0。

独立证据适配器 v1 把原生 `sha256:<hex>` 与裸十六进制哈希直接比较，失败关闭；v2 在 QC review 的线粒体比例数值核验处失败关闭；v3 在 NumPy 字符串与 JSON 字符串的严格类型比较处失败关闭。交并、Jaccard、匹配标签实际值已交叉核对一致；身份字符串显式归一后，v4 全部通过。三份失败收据、日志及对应源码快照均保留，原生制品、包源码、冻结协议与运行环境没有因此修改。

数值诊断及 v4 同时要求整数归约/比例分子精确、源 dtype 公式与 CSV 差不超过 1e-12、数学 float64 值与最终 float32 除法的差不超过半 ULP 加 CSV 误差界。PBMC 比例最大数学差为 6.640229e-7 个百分点，源 dtype 路线最大差约 1.78e-15；完整 14 条运行的缺失位置和阈值判定一致，没有阈值翻转。Baron 无 MT 基因，比例正确保持 NA。浮点舍入界不用于放宽整数计数、轴、名单、分区或指标比较；ARI/NMI 等指标重算容差仍为 1e-12。

原生总耗时的 monotonic 起点仅保存在失败 launcher 的内存里，没有有效持久收据，因此本轮不报告原生总耗时、峰值内存或性能提升。CI 自身的有效耗时收据保留；它们不是原生科研流程的资源测量。

## 下一轮应改进的具体环节

1. **把长运行与最终汇总分开提交状态。** 在每个阶段开始/结束时持久保存时间与完成状态，最终 formatter 失败应明确标记 postprocessing failure，并从既有原生结果恢复；本轮外层失败证明内存里的计时信息不足以恢复资源证据。
2. **为证据适配层建立类型约定。** 统一有算法前缀的 digest 与文件裸哈希的边界，明确统计量的源 dtype、计算路线、单位、缺失规则与可证明舍入界，在大型冻结运行前用实际 H5AD/CSV fixture 验证适配器。数学整数计数、细胞/特征轴、名单和阈值判定仍应精确一致。
3. **按局部群体与供体组织后续科学问题。** 总体 ARI、簇匹配、完整主簇覆盖和外部标签一致性分别返回。对留出后仅剩少量细胞的群体明确不可据此估计跨供体稳定性；需要新的供体/数据与提前固定的成功准则来支撑确认性重现，不能通过同数据选最高分分支代替。
