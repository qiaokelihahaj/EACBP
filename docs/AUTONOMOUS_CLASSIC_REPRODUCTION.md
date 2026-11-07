# 经典单细胞项目的自主复现与架构改进

本轮遵循用户要求：公开实验数据、下载依赖、缓存、临时文件和分析结果放在 `F:/EACBP-autonomous-practice`；EACBP 自主接收科学问题与原始输入，选择已实现的方法和分析分支。人只选择教学项目、提供物种/组织背景和计算资源，独立评价器在原生分析完成后读取参考答案。

## 自主执行的范围

入口为 `eacbp analyze --question ... --data ... --species ... --tissue ... --run-dir ... --threads 2`。调用者不传 QC 阈值、目标细胞群、预期 marker、轨迹根细胞或科学 workflow 配置。该入口是依据已观测矩阵和实验设计的确定性规则系统，尚不是通用语言模型科学推理服务。

项目识别原始整数 counts、条件和 donor 别名，检测配对/样本资格，依据观测分布生成稳健 QC 阈值。探索问题进行原生 Scanpy 预处理、聚类和描述性群体特征排名；明确要求比较已观测群体的条件问题，对所有有资格群体执行 donor pseudobulk 和完整 leave-one-donor-out；请求 RNA velocity 且存在合法 spliced/unspliced 层时调用真实 scVelo。缺少必需方法、未估计模型或缺失 donor fold 不能通过问题完成验收。
HVG
执行使用安装后的冻结镜像、非 root 用户、只读文件系统和禁用的网络。只挂载当前原始 H5AD、输出、依赖和缓存；参考目录没有挂载。执行问题与完整命令保存在每次运行的 `*-execution-command.json`。所有失败运行和原始日志保留，重跑创建新目录。

## 固定输入与独立评价
维度 | SHA-256 |既有
导入时保留 counts、基因身份和实验设计。历史聚类、嵌入、邻居图、 cell_type 是问题明确要求使用的观测分组，并不代表系统自主识别了这些类型。
3ac27a2a2d12581fcd77afed7ea653a1` |
| Kang 2018 IFN 响应 | `kang_2018.h5ad` | 24,673 × 15,706 | `e6a5adac64dcdeb36eaba27db49b63e0c64bb0ed4a64c6705971506b41c39830` |
| 小鼠胚胎胰腺 RNA velocity | `endocrinogenesis_day15.h5ad` | 3,696 × 27,998 | `9e3e459eca00ba06b496ec80def32941b5b2889918720e3e7aa6ffb811fbe7
|---|---|---:|---|
| PBMC3k 群体探索
| 教学项目 | 公开输入 |  | `pbmc3k_raw.h5ad` | 2,700 × 32,738 | `89a96f1beaa2dd83a687666d3f19a451、根/终点及其他派生答案被隔离；胰腺的原始剪接层保留。Kang 的c6` |

来源说明和接收记录保存在 `references/source_receipt.json`。胰腺输入来自 [scVelo 官方数据接口](https://scvelo.readthedocs.io/en/stable/scvelo.datasets.pancreas.html)及 [Velocity Basics 教程](https://scvelo.readthedocs.io/en/stable/VelocityBasics.html)，关联 Bastidas-Ponce 等的 [2019 年研究](https://doi.org/10.1242/dev.173849)。该胰腺教学 H5AD 没有单独声明的数据许可，不能将软件许可当作数据许可。接收记录另外保留 PBMC 上游 10x 数据及 Kang Figshare v2 的 CC BY 4.0 声明。

独立评价方案 `references/blind_acceptance_protocol.json` 的 SHA-256 为 `184a912820a0a8e165cc642b4bb2f1484727f0af3d920a967ee7d12409502d02`。方案在本轮原生新结果产生前冻结；此前已运行的手工配置验证属于已知背景，因此这不构成对全部历史的事前注册。科学阈值没有根据本轮结果调整。

评价器先检查原生问题完成状态、运行/制品/审计哈希及执行隔离，再读外部参考；不改原生配置或产物。工程完成与科学复现通过分别报告：

- PBMC：总体与逐类保留比例、历史群体的 ARI/NMI、六类描述性 marker 覆盖；总体得分不能覆盖稀有群体失败。
- Kang：所有合格观测群体、完整 donor 配对、IFN 基因方向、独立 donor 表达方向、关键群体 FDR 与 donor 一致性、全部 LOO 拟合完成情况。巨核细胞因 donor 集合不完整而不合格，资格来自原始设计，不来自参考结果。
- 胰腺：总体与逐类保留、群体结构、速度基因与模型置信度，以及 Ngn3 high EP → Pre-endocrine → Beta 的全局速度图方向。方向评价用按源细胞归一化的正/反向边均值，要求足够跨状态边；分析阶段不传这些标签或根节点。

评价器的追加复核修复了几个边界漏洞：首次读取参考前验证完整 snapshot、summary 与任务/审计；执行记录绑定本次目录、原始文件实际哈希、问题/背景和独立核对的镜像摘要；解析挂载的绝对路径、祖先及 junction/symlink，拒绝未知形状；速度 CSR 合并重复并去零，仅计唯一正权跨群边。追加 **34 项测试全部通过**，证据为 `runs/evaluator-hardened-final.*`。本地执行 receipt 的一致性核验不是签名的远程执行证明。

## 原生运行与复现结果

PBMC 的第一轮原生运行 `runs/pbmc3k-autonomous-v1` 已完成全部六项必需任务并通过审计。独立评价 `runs/pbmc3k-autonomous-v1-evaluation.json` 报告工程通过、完整科学复现未通过。与历史 processed reference 共有 2,557/2,638 个细胞（96.93%），ARI 为 0.5633、NMI 为 0.7388；六类主要群体的特征基因达到固定标准，但巨核细胞仅保留 5/15（33.33%），低于逐类 80% 的标准。DC 虽全部保留，也没有形成高纯度独立簇。系统选择的 400 基因 QC 下限没有改为参考教程的手工阈值。

Kang 的 v1 运行保留为失败证据：24,673 个细胞经过合法 QC 保留 24,555 个，保留细胞的标签没有改变；保护器将删行误判为重分群。胰腺的 v1 运行也保留：扩展工厂错误地将尚未经过实际方法路由的模板方法作为最终方法检查，执行前就拒绝任务。两者的 v2 使用冻结修复镜像 `sha256:ab952586373f71eca2076c8b6b9ecdd4f87665a86591d4dc6e340db2c51128dc`，源文件核对记录为 `runs/installed-source-v2-realfix.json`。

胰腺 v2 已完成七项原生任务，真实 stochastic RNA velocity 与独立审计均通过。独立评价 `runs/pancreas-autonomous-v2-evaluation.json` 沿用冻结方案，报告工程与完整科学复现均通过：保留 3,695/3,696 个细胞，各类保留至少 99.84%；ARI 为 0.4454、NMI 为 0.6804；1,658 个速度基因，速度置信度中位数 0.8703。Ngn3 high EP → Pre-endocrine 和 Pre-endocrine → Beta 的方向指数分别为 0.5546 和 0.4665，对应 3,648 和 7,880 条跨状态边。分析没有接收参考群体或根节点。这证明固定教学范围内的方向一致，不确认真实命运或因果。

Kang v2 的 27 项计算和原审计均已完成，覆盖全部七个合格群体、七个主模型、七个 LOO 基线重拟合及 56 个 donor 删去拟合。随后展示保存因 CSV 空描述被解析为 NaN 而失败；原事件日志仍以 27 项完成、零阻断的 success 结束。恢复只从该已完成 journal 重建报告，不重拟合、不重新出具科学审计。

实际恢复使用镜像 `sha256:c183a4930ed9aa6bcf4ca1a30fa378e7e608709f01f739684434f045af364279`，安装源码核对记录为 `runs/installed-source-v5-admitted-warning-recovery.json`（118 个源码文件一致）。报告成功发布，问题状态为 complete；原计算镜像和环境保持 v2。恢复回执 `runs/kang2018-autonomous-v2/finalizations/9d543f0232054b42ab6e5ad1c223caa6/receipt.json` 的 SHA-256 为 `8fe8151b616e8e1bb8ee54f058ad7f33922bd321175259e3245c88ca49303489`，记录零科学执行、零新科学审计及原/现版本。

外部评分前的独立核验 `runs/kang2018-native-coverage-verification.json` 已通过，SHA-256 为 `422505ddaea2d340b4ce790d8f9ccd1ea3239c84c15c6f5840fe62c0a7040c02`。原 100 个文件除允许发布的 summary 外均未改变，失败摘要的逐字节备份一致；所有科学产物、原签名、审计和日志保持。每组 8 个配对 donor、16 个样本；全部 56 个 LOO fold 具有完整基因集合，删去对应的两份样本后为 14 个样本、每条件 7 个 donor。保留细胞的原始 donor/condition/population 标签一致。核验器没有读取参考答案；恢复回执中的零调用声明由冻结恢复代码、禁止计算调用的测试及原科学记录不变共同支持，不是签名的远程执行证明。

Kang 的最终独立评分 `runs/kang2018-autonomous-v2-evaluation-hardened.json` 报告工程通过、完整科学复现失败，唯一失败项仍是逐类保留不足。总体保留 24,555/24,673（99.52%），但巨核细胞仅保留 22/132（16.67%），低于该项目固定的 70% 逐类门槛；这项保留要求涵盖所有原始观测群体，不因该群体缺少完整配对 donor 而豁免。其余七个群体的预定七个 IFN 基因全部方向为正、FDR < 0.05，独立 donor 表达方向也全部为正；两个预定锚点 B 和 CD14+ 单核细胞的 donor 一致性及 FDR 标准通过。全部合格目标和 LOO 覆盖通过。评价前后原生文件哈希一致，方案和科学参数没有改变。

| 项目 | 原生问题完成与工程验收 | 固定范围内已复现的结果 | 完整外部科学验收 |
|---|---|---|---|
| PBMC3k | 通过 | 主要群体结构与六类特征基因 | 失败：巨核细胞保留 33.33% |
| Kang IFN | 通过，展示失败已安全恢复 | 七个合格群体的 IFN 响应及完整 donor/LOO 分析 | 失败：巨核细胞保留 16.67% |
| 胚胎胰腺 | 通过 | 群体结构及两段预定 RNA velocity 方向 | 通过 |

新增 QC 可观测性的 PBMC 复跑为 `runs/pbmc3k-autonomous-v2-qc-records`，使用冻结镜像 `sha256:7b3d299dbd3910ca24418f8472e63b7893a2f5187665103fc26ee9b3622f5c34`。新 QC 产物覆盖全部 2,700 个输入细胞，独立审计与只读重验通过：保留 2,562 个，仅低基因数 81 个，仅高线粒体比例 41 个，两项同时触发 16 个。与 v1 的六个任务科学结果逐项完全一致，包括 counts、PCA、UMAP、邻居图、聚类和 180 × 17 marker 表，所有数值差为零；阈值未改。核验为 `runs/pbmc3k-qc-observability-verification.json`，独立科学评分 `runs/pbmc3k-autonomous-v2-qc-records-evaluation.json` 仍仅因逐类保留不足而失败。

最终以加固后的评价器重新评分，PBMC 与胰腺的科学指标和结论均未改变。新记录为 `runs/pbmc3k-autonomous-v2-qc-records-evaluation-hardened.json` 与 `runs/pancreas-autonomous-v2-evaluation-hardened.json`；此前评价保留供追踪。

包含报告恢复入口、CSV 描述语义修复及非阻断审计提示比对的最终源码，核心 CI 为 **565 passed、1 明确 deselected、无未预期 skip**；standard lane 为 **154 passed、1 明确 deselected、无未预期 skip**。JUnit 与日志为 `runs/core-final-v5-admitted-warning-recovery.*`、`runs/standard-final-presentation-recovery.*`，门禁已实测通过。standard lane 在最后仅影响独立恢复入口的提示比对修复前运行；该入口的最终专项测试为 **35 passed**，已纳入最后核心 CI。Linux NumPy 2.3.5 环境此前的真实 scVelo 与相关回归为 **40 passed**，证据为 `runs/velocity-compatible.*`；真实胰腺原生运行进一步覆盖本轮方法路由修复。不同检查可能重叠，不将计数相加为独立测试数。

## 已落实的架构改进

| 实验暴露的问题 | 通用修复 | 适用边界 |
|---|---|---|
| 流水线 success 与科学问题完成混淆 | 持久化问题契约并核验必需任务、目标分支、可估计模型、审计及所有 donor fold | 任务完成仍不代表外部科学复现通过 |
| 旧 H5AD 已带结果可能污染自主分析 | 导入重新生成原始 counts/设计输入，隔离历史派生答案，并保存准备 receipt | donor、条件和既有观测分组的实验含义仍需真实元数据支持 |
| 探索流程没有可审查的群体特征 | 注册描述性 cluster marker 表，审计独立重算排名/数值/覆盖 | 同一批细胞定义聚类与 marker，不是独立假设检验或自动类型认定 |
| RNA velocity 缺少原生能力及审计 | 注册真实 scVelo stochastic 分支、输出摘要、counts/邻居/moments/数值一致性审计 | 不独立重拟合全部动力学，不确认生物学假设；方向是模型推断 |
| NumPy 2.4 以上破坏 scVelo 0.3.4 标量转换 | 限定兼容版本并锁定实际 NumPy 2.3.5；不使用模拟速度或库补丁 | 新库版本需要另行兼容验证 |
| 可覆盖的依赖列表错误记录基础环境版本 | 通过活动 metadata 搜索路径解析各 distribution 的版本，并与本轮实际导入版本核对 | 本轮安装/执行版本已核对；该 resolver 本身不独立证明任意导入路径或跨版本按位一致 |
| 方法别名导致有效审计回执不能用于完成判断 | 通过实际能力注册表规范化方法身份，并继续校验契约/输出指纹 | 未注册方法、错误方法及篡改指纹仍被拒绝 |
| 合法过滤被误判为标签改写 | 按保留且唯一的 barcode 比较标签含义，兼容 H5AD category roundtrip | 真正标签改写、缺失标签、未知/重复 barcode 仍被拒绝 |
| 扩展规划误读尚未路由的模板方法 | 工厂获得按实际 router 解析的深拷贝上下文；原始契约保持不变 | 真正非原生 Scanpy 聚类仍不能作为 RNA velocity 邻居来源 |
| 全局 QC 可能丢失低 RNA 稀有群体 | 原生 QC 增加全输入细胞处置记录、独立审计和解释限制 | 阈值/保留 mask 不改；自动保留低 RNA 的敏感性分析分支尚未实现 |
| 已审计计算完成后，CSV 空描述被解析为 NaN 导致报告保存失败 | 按已审计元数据恢复空 covariates/status_reason 的语义，证据及快照序列化仍严格拒绝数值 NaN/Inf | 原统计表允许有明确未估计状态的 NaN；不清洗统计缺失值、不改原始结果表 |
| 展示阶段失败使已完成的昂贵计算难以恢复 | 独立的已完成 journal 报告恢复入口，核对原环境/签名、完整任务、制品和审计，再重建展示 | 普通 resume 的源码约束保持；原计算与当前展示版本分别记录 |
| 事件里的非阻断诊断被误当作审计拒绝 | 事件的失败检查名单必须与原持久审计报告完全一致，保留 warning/info，继续拒绝 error/stop_rule 和未通过回执 | 提示信息不会被删除，也不会使不合格科学产物获得准入 |

稀有群体失败说明全局分布阈值不能独立判断所有低 RNA 细胞。参考 [OSCA 质量控制讨论](https://www.bioconductor.org/books/release/OSCA.advanced/quality-control-redux.html)，本轮增加可核查的处置记录，但不按已知参考类型反向调整阈值。未来需要依据独立实验信号进行自动 QC 敏感性分析，再用新的数据评价。

NumPy 的兼容边界还由 [2.4 官方发行说明](https://numpy.org/doc/2.4/release/2.4.0-notes.html)确认：非零维数组到标量的旧转换开始直接报错。版本限制来自实际失败和上游 API 变化，不来自速度方向的预期结果。

## 复跑

在仓库根目录运行，使用新输出名称，不覆盖历史结果。`--run-name` 仅改变保存目录，不选择科学参数：

```powershell
.\.venv\Scripts\python.exe -B scripts/run_classic_autonomy_curriculum.py `
  --root F:/EACBP-autonomous-practice --image <已核对的冻结镜像 ID> `
  --project kang2018 --run-name kang2018-new-replay --threads 2

.\.venv\Scripts\python.exe -B scripts/evaluate_classic_autonomy.py `
  --project kang2018 --run-root F:/EACBP-autonomous-practice/runs/kang2018-new-replay `
  --expected-image <同一已核对的冻结镜像 ID> `
  --reference-root F:/EACBP-autonomous-practice/references `
  --output F:/EACBP-autonomous-practice/runs/kang2018-new-replay-evaluation.json
```

对全部科学计算与原审计已完成、但展示保存失败的目录，使用独立恢复入口。它拒绝不完整、篡改或已经成功发布的目录；普通 resume 不承担跨源码版本的展示恢复。Kang 本次实际执行的完整镜像、挂载和参数保存在 `runs/kang2018-finalization-command-v5.json`。其原生命令为：

```text
eacbp finalize-completed-journal
  --run-dir /results/kang2018-autonomous-v2
  --compute-image sha256:ab952586373f71eca2076c8b6b9ecdd4f87665a86591d4dc6e340db2c51128dc
  --source-proof /proof/installed-source-v2-realfix.json
  --source-proof-sha256 62b0b46934fa998f57dda35b5fcc630378853dc5f03ecde76f1d9e1b6b601f12
```

科学结论仅限于这些公开教学项目。Kang 原研究的 demuxlet/遗传分析没有复现；速度和 pseudotime 不是观测的真实时间、命运或因果；同一数据内审计、marker 与 LOO 不能替代独立队列验证。未估计、失败及保留不足都作为结果记录。

## 本轮保留的局限

除了 PBMC 稀有群体的 QC 敏感性，真实运行还暴露两个可继续改进的工程观察。胰腺全基因多层 RNA 产物约 2.087 GB，执行/保存约 211 秒、提交约 51 秒、独立审计约 243 秒；同一不可变大产物在多个检查入口重复加载/哈希。后续可通过同次审计内共享只读载荷和已验证摘要减少开销，当前没有减少科学检查或偷偷下采样。

Kang 的运行日志出现 PyDESeq2 参数化 dispersion trend 不收敛、自动切换 mean-based trend 的警告。这是库的实际拟合行为，不能描述为全部参数拟合均收敛。保留的合并日志没有每条警告的任务/折标识或时间戳，receipt 尚未记录逐拟合 effective trend/warning，因此不能据警告行数推算受影响任务、donor fold 或失败拟合数。这项可观测性缺口保留为后续改进，没有按预期结果修改模型设置。

冻结的速度图指标比较按源细胞平均的跨群正向/反向流量，受群体大小和邻域结构影响，不是无偏的成对方向统计检验。本次核心边的群体大小分别为 641/592 和 592/591，仍以原先定义和阈值评分；通过只表示固定跨群流量标准与教程方向一致，不替代基因动力学相图审查或实验谱系验证。
