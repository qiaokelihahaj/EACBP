# P4：长任务恢复、阶段收据与新公开数据验证

本轮实验数据、下载、缓存、临时文件、日志和证据位于 `F:/EACBP-architecture-p4`。源码与本文件保留在工作仓库。P3 索引的 135 份原绑定源码已逐字节复制到 F 盘历史源码目录，P1/P2/P3 的旧索引、文档和原生结果保持原状；P4 修改当前工作区源码，不宣称旧索引指向的当前源码仍与 P3 相同。

## 已实现的恢复边界

`execute_clustering_sensitivity(prepared, resume=True)` 显式恢复预先冻结的实验。身份绑定源文件、plan、bindings 和实际包/依赖环境。每次 invocation 与 preflight、native_execution/native_admission、comparison、postprocessing 分别保留开始和结束事件。事件以完整文档独占发布并 fsync；head 为原子检查点，能发现尾事件丢失，并从“事件已完整发布、head 尚未更新”的崩溃点恢复。进程锁防止同一实验并发执行。这是本机一致性检查，不是外部签名证明。

恢复时先核对历史已接纳分支的文件名单/哈希，再重新调用原生 collector，检查全部持久审计、制品、计数、参数、manifest 与轴顺序。状态为 completed 本身不能授予复用。已接纳目录丢失或字节变化会阻止恢复，不能把它当新分支重新计算。没有开始的分支可以继续；已有但不完整的原生研究保持原状，任务级恢复仍使用 `study_service.resume_study`，此 API 不自动重试它。

崩溃可能发生在原生完成之后、首次 collector 之前。此时 sensitivity 日志没有旧原生文件 pin，恢复必须依靠真实持久原生任务证据与冻结输入重新接纳，并明确返回 `resumed_without_previous_native_file_pin`。本轮故障实验另在进程退出前独立保存原生文件哈希，以验证该实际案例的字节保护。

新 invocation 生成自己的 comparison/report 版本，旧版本保留；`latest_comparison.json` 指向最近成功发表的版本。首个根目录 comparison/report 作为历史兼容入口保留。恢复后的调用者应使用 API 返回路径或 latest 指针，不能把首个历史文件当最新结果。恢复 API 强制重核历史原生证据；旧 detached comparison/report 不是当前科研输入，也不重新准入，其外部改写不在 API 的恢复检测范围内。本次故障实验另独立核验旧比较版本在恢复前后字节不变。

报告失败记录为 postprocessing_failed，已经完成的原生接纳与比较阶段保持完成记录。正常同进程阶段有持久 monotonic 耗时；旧进程未结束的尝试追加 interrupted，耗时为 null，不能用墙钟差补造计算耗时。native_execution 的计时包括运行与接纳，native_admission 为恢复复核；两者不等于纯算法 CPU 时间。

## 资源观测

每个 invocation 在实际 Python 进程中记录进程树采样，逐行 fsync JSONL，正常结束另写原子摘要。Windows 使用 Toolhelp32、进程创建身份和 WorkingSetSize；Linux 实现 `/proc` 路径，本轮本机实测范围为 Windows。记录已观测的真实科学子/孙进程，而非仅测 venv 启动器。PID 重用与探测错误单独记录。

`observed_peak_sum_rss_bytes` 是顺序、非原子地读取各进程 RSS 后所得 SUM 的最大观测值。共享页可重复计数，瞬时峰和未观测到的孤儿进程可遗漏，因此它不是物理内存峰值、数学保证的下界、完整 job 资源核算或生产内存上限。实际最大采样间隔另列；设定 interval 不等于准确周期。没有成功观测时值为 null，coverage 为 unavailable；有错误时保留 partial。摘要格式化失败时已有采样日志保留，不掩盖原科学异常。

## 新数据与预先声明的边界

Muraro GSE85241 的完整公开表含 UMI 校正后非整数值，严格整数计数入口拒绝它；原文件、判定与日志保留，没有舍入后冒充原始计数。

新候选采用 Segerstolpe E-MTAB-5061 官方 deposited processed 文件的 raw-read-count 部分，完整保留 3514 个细胞及 26271 个原始 feature 行。IDF 的 26272 是含 header 的物理行数；首次提取因错误解读该行数而拒绝输出，失败日志和源码保留，没有截断源数据。独立 C parser 逐值核对全部 92,316,294 个 read counts，零容差相等。RPKM 部分不进入计数矩阵。重复 gene symbol 通过唯一 feature 标识保留全部行；外部 marker 描述预先规定按原 symbol 对所有对应行求和。93 个 ERCC/eGFP feature 也保留，完整源轴与生物基因范围不是同一个概念。原 symbol 的 MT- 前缀匹配为 0，线粒体指标为 unavailable，不能记作 0 或解释为没有线粒体转录本。

只有真实 SDRF cell ID 与 individual 字段参与原生输入，10 名 donor；作者 cell type、inferred cell type 与质量列保持在独立参考文件中。H1 有 96 个细胞，且官方说明其 ERCC 加样量为其余样本四倍；这是已知技术混杂，当前默认流程保留全部 feature，不能把原始来源完整性当成技术混杂已经纠正。数据包含健康与 T2D 供体，元数据不建立随机抽样独立性或临床因果结论。

此前 Baron resolution=0.5 是探索性发现。本轮在新数据分数读取前固定 resolution=1 主分析与 0.5、2、seed=17、原始细胞 80% 描述性分支。10 donor 超出既有 2–4 donor 全 LOO 策略，LOO 明确不可用，不能声称已经验证跨供体稳定性。

新数据原生 QC 沿用已有完整输入上的 3 MAD 规则，冻结为 min_genes=689、max_mito_pct=100；MT 标识不可用时 100 是占位阈值，不能作有效线粒体筛选证据。准备草案曾提出固定 400/20，但首次冻结因遗漏当前源观测绑定而被拒绝；失败日志和草案源码保留。修正准备脚本、恢复原有策略并冻结后才开始原生分析，未读取新分数作阈值调整，包内科研策略没有修改。主参考假设为 resolution=0.5 相对 resolution=1 的 ARI 增加至少 0.05，同时 ARI≥0.5、NMI≥0.6、完整有标签源细胞整体及每类保留率均≥0.8。分数、覆盖和假设检验分别返回；不可估计不记作成功。阈值不会因实际结果改写。归一化、2000 HVG、30 维上游表示、15 neighbors 和整型源计数范围沿用已有流程；ERCC 与 eGFP 对默认流程的影响仍属于限制。

当前官方 deposit 的 inferred cell type 预先指定为主参考字段，cell type 为次参考，两套分别报告，不能看结果后择优。单类、没有足够共同细胞或没有有效标签时相关分区指标不可估计，返回 null。作者分类依赖表达数据，当前 deposit 与原论文 final figure 标签的版本关系没有独立证明；结果限于明确来源字段的一致性，不等于原论文算法、图或独立生物真值的完整重现。

[原研究实验室的说明](https://sandberglab.se/pancreas/)报告筛选后的 2209 个细胞及七种主要细胞类型。本轮读入完整 deposited 3514 个细胞，未用作者质量标签挑出 2209 细胞子集，两者分析总体不同。完整源保留便于查清排除分母，不授予原论文图结果已经复现的结论。

获取下载地址时，acquisition 脚本曾整体解析并保存含作者 cell type、inferred cell type、quality 属性的官方 BioStudies API JSON。因此不能声称全团队从未解析标签或完全盲确认。SDRF 分类值在原生冻结、执行与独立接纳完成之前保持关闭，冻结与验证补充声明创建时新的外部分数均未读取或计算；这是披露元数据暴露的预先规定检验。公开元数据、技术路线与此前 Baron 评分属于已知背景，应与新的参考内容何时读取分别记录。

## 完整回归与已完成的恢复实验

完整 CI 共 1267 项通过，0 failure/error/skip：core 959、standard 237、advanced-statistics 71。保留此前 core 与 standard 各一项明确 deselect，分别为跨进程 H5AD 全 slots 回载和 DPT→CellRank 计划连接。没有增加跳过项或为本轮扩大排除。

新增单独测试包含 journal 23、进程树观测 16、恢复集成 22；独立数学 helper 24、native adapter 11 类、事件/资源/故障证据 27、外部科学 checker 51 项合成案例均通过。它们与 CI 或定向回归可能重叠，不累加成独立科研验证数量。

| PBMC 阶段 | 实际子进程退出码 | 新原生运行 | 已有原生重新接纳 | 子进程墙钟秒 |
| :--- | ---: | ---: | ---: | ---: |
| primary 返回后强制退出，首次 collector 前 | 23 | 1 | 0 | 59.3117 |
| 首次显式恢复，最终 formatter 失败 | 1 | 4 | 1 | 185.0239 |
| 最终显式恢复，禁止 run_study/resume_study | 0 | 0 | 5 | 24.8075 |

这里的墙钟来自外层完整子进程观测，包含 Python 启动、输入核验、读取、写出及原生运行，不能当成纯算法耗时或中断阶段的恢复耗时。首次进程的 invocation/native_execution 仍由后续恢复记为 interrupted、elapsed=null；有采样前缀但没有完成的资源摘要。

PBMC 本轮 primary 与 P3 primary 的实际持久准入输出已作零容差比较：QC、归一化、integration 表示、聚类、marker 表的科研数组、表格与轴顺序相同。URI、study ID、时间与诊断元数据不作为科研数值等价；没有用字节不同来要求重算旧科研结果。

真实 prepared binding 的独立只读校验已通过全部 10 个范围：PBMC 完整源 2700 细胞，固定 QC 保留 2562；Segerstolpe 完整源 3514 细胞，冻结 3 MAD 基因数阈值保留 2933。QC 候选与计划抽样遗漏分别记录，未把候选恢复到主分析。

## 数值存储核验发现与修正范围

首次完整独立 native verifier 在 Segerstolpe primary 的描述性 QC feature mean 表失败并关闭外部评分入口。105084 行中 52542 个有限均值有两处经 pandas 默认 CSV parser 解码后超过原有绝对界 1e-12，最大差 1.8189894035458565e-12。例如 CSV token `12538.393112853733` 被默认解码为 `12538.393112853732`。较小均值也可能相差较多 ULP，不能用全局“只允许 1 ULP”掩盖真实解码范围。

独立只读诊断证明：原始 CSV token 经 Python float 和 pandas round_trip 两条路径解码，全部 52542 个有限均值与精确整数归约后源 dtype 除法的结果零容差一致；实际 get_audited 返回的默认解码表也与默认 parser 重放逐值一致。QC disposition、计数、轴与 retained mask 均未改变，描述均值不参与 QC 阈值判断。

修正保存为新的 `verify_native_stability_v2.py`，先核实持久准入与默认解码载荷，再对字节绑定 CSV 的原始数值作精确公式证明，将默认解码漂移单独返回。计数和轴仍零容差、分区指标仍 1e-12、原有 float32 最终舍入证明保持原状。旧 v1 checker、失败收据、协议、mask、supplement、包源码、原生结果都未修改，原生重跑次数为 0。

这是原生执行后、外部评分前的验证补充声明，不冒充第一次原生冻结前的 checker。评分入口另保存为 v2，使用实际 v2 身份；25 项科学 AST 投影，包括参考解析、科学常量、完整评分循环和主假设赋值，与旧版一致。71 项测试中复用原 51 项科学案例，另加 20 项身份、补充声明和 AST 漂移拒绝案例；CSV 证明另有 6 类聚焦验证。原版失败与新版通过应同时保留，不能只返回成功摘要。

## 实际资源收据

| invocation | 采样条数 | 最大观测 RSS SUM 字节 | 实际最大采样间隔秒 | 完整子进程墙钟秒 |
| :--- | ---: | ---: | ---: | ---: |
| PBMC 强制退出 | 148 | 834211840，仅持久前缀 | 2.7644 | 59.3117 |
| PBMC 恢复后报告失败 | 459 | 1064169472 | 2.7354 | 185.0239 |
| PBMC 最终恢复 | 59 | 706850816 | 0.4079 | 24.8075 |
| Segerstolpe 五分支原生执行 | 2122 | 4690702336 | 2.9874 | 775.9480 |

设定采样周期为 0.25 秒，实际间隔明显更大，必须随结果返回。四次原生 invocation 观测到的进程最大数均为 1；这不证明进程树已完全覆盖。进程树功能另外用真实分配内存的科学孙进程测试验证。中断前缀没有最终摘要，也没有完整 invocation monotonic elapsed；另三次属于 sampled_observations。约 4.69 GB 的 SUM 是本次样本统计，不能据此承诺生产内存限制或推断未测试平台性能。

P3 的 851 份旧证据、135 份原源码归档及其他旧绑定记录合计去重 1009 个保护文件，在本轮全部原生执行之后逐字节核验无变化。当前工作区源码按 P4 演进，旧版本通过精确归档追溯。

## 完整独立稳定性结果

v2 通过全部 10 个原生范围、70 次任务持久准入及 8 个配对的独立重算，完整源计数、source dtype、轴、donor、原始抽样 mask、QC 分区和上游 embedding/图/Leiden 收据一致。原生成功与 detached 科学评分分开，不把外部假设是否通过写回 completed 状态。

| 数据 | 分支 | 共同保留细胞 | ARI | NMI | 一对一未匹配主簇数 |
| :--- | :--- | ---: | ---: | ---: | ---: |
| PBMC | resolution 0.5 | 2562 | 0.837731 | 0.878536 | 3/10 |
| PBMC | resolution 2 | 2562 | 0.703906 | 0.841896 | 0/10 |
| PBMC | seed 17 | 2562 | 0.927130 | 0.926658 | 0/10 |
| PBMC | 原始源抽样 80% | 2046 | 0.844114 | 0.878746 | 1/10 |
| Segerstolpe | resolution 0.5 | 2933 | 0.845905 | 0.931076 | 7/30 |
| Segerstolpe | resolution 2 | 2933 | 0.782448 | 0.932603 | 0/30 |
| Segerstolpe | seed 17 | 2933 | 0.924749 | 0.958178 | 0/30 |
| Segerstolpe | 原始源抽样 80% | 2358 | 0.878502 | 0.931268 | 2/30 |

未匹配表示最大总 Jaccard 一对一分区匹配中没有对应簇，不等于细胞消失或已证实的生物群丢失。全局高一致性仍可能掩盖局部变化：PBMC seed 17 的 ARI 为 0.927130，但最小主簇 Jaccard 为 0.147541。80% 源抽样后完整主簇的最小共同覆盖分别为 0.756701 和 0.714286；分母是完整主簇，不能换成抽样子集以放大覆盖。每个主簇与未匹配簇均保留在结果中，小支持少于 20 的匹配单独标记。

Segerstolpe donor 覆盖按真实来源、branch raw 和 retained 三层独立重算；10 donor 的 LOO 仍不在本轮范围。这八对是同数据相关敏感性分析，不能作为八次独立生物复制。

## 新参考字段的实际科学结果

原生 v2 全部接纳与 8 对独立重算通过之后，评分入口重新准入全部 10 个范围，才读取真实 SDRF 的两套分类值。主字段有 2209 个有效标签细胞、14 个原始类别，另 1305 个细胞标签缺失；源 ID 与 SDRF 的 3514 条记录完整一一对应。有效类别仅 strip，不按大小写合并、不把 unknown/unassigned 删掉或改作 missing。次字段完整源全部 3514 细胞只有一个类别，五个分支 ARI/NMI 均不可估计并返回 null。

| Segerstolpe 分支 | 主参考共同有效细胞 | 主参考 ARI | 主参考 NMI | 有效标签完整源保留率 | 主参考验收 | 次参考分区指标 |
| :--- | ---: | ---: | ---: | ---: | :--- | :--- |
| primary，resolution 1 | 2209 | 0.277308 | 0.656911 | 1.000000 | 失败：ARI | 不可估计 |
| resolution 0.5 | 2209 | 0.327864 | 0.683832 | 1.000000 | 失败：ARI | 不可估计 |
| resolution 2 | 2209 | 0.238612 | 0.644568 | 1.000000 | 失败：ARI | 不可估计 |
| seed 17 | 2209 | 0.286491 | 0.670142 | 1.000000 | 失败：ARI | 不可估计 |
| 原始源抽样 80% | 1787 | 0.303304 | 0.684433 | 0.808963 | 失败：ARI、每类保留率 | 不可估计 |

主假设的 ARI 增量为 0.05055599759377144，达到 0.05；但 target ARI 为 0.3278644226787174，未达到 0.5，因此组合假设明确为 failed。NMI、整体有效标签保留率与每类保留率通过不能替代 ARI 条件。五个主参考分支均未达到完整绝对验收，次参考不可估计也不算通过，没有从两套字段择优或选取新默认。

前四个分支都保留全部 2209 个主字段有效标签细胞，没有计划源遗漏、也没有这些细胞的 QC 候选。标签缺失的 1305 细胞中保留 724、QC 候选 581，缺失不被解释为低质量真值。80% 分支对有效标签细胞的 422 个遗漏全部来自冻结源抽样，QC 候选为 0；缺失标签细胞另有 281 个抽样遗漏、453 个 QC 候选与 571 个保留。两种排除原因不合并。每个原始参考类别，包括小类别与低覆盖类别，都留在 JSON 的完整源分母中。

有效标签数 2209 与论文实验室说明的筛选细胞数相同，但本轮没有逐条证明当前 deposit 的分类版本与原论文 figure 标签完全对应；14 类也不等于七个主要类型。结果证明的是在固定来源字段上的一致性及冻结假设失败，不授予完整原论文图复现、独立生物真值或新发现的结论。

80% 分支每类保留率失败的原始字段类别为：MHC class II cell 2/5、alpha cell 697/886、endothelial cell 12/16、epsilon cell 5/7、gamma cell 154/197、mast cell 5/7。unclassified cell 的 2/2 与 unclassified endocrine cell 的 34/41 也保留在同一分母中，没有事后合并稀有类或删除未分类类。

## 后续改进应由这些证据驱动

数值制品应在下一轮引入明确的数值类型和存储/解码合同，并用实际大型 CSV 保存与持久准入载荷作为回归 fixture。此次补充证明解决验证边界，没有修改当前包的默认 CSV 读写；未来更改需要新的源码冻结和科研数值回归。

恢复能力已有应用层 API，下一步可以把显式恢复及最新结果路径接到 CLI/WebUI，并返回可直接阅读的“已有任务重新接纳、未开始任务继续、失败阶段保留”清单。不完整原生研究的任务级恢复与旧报告篡改检测仍应分别定义，避免一个 resume 标签承诺过多行为。

结果呈现应同时保留工程状态、科学假设状态、参考来源、完整分母、不可估计项、每簇局部变化、失败版本和修正声明。证据索引必须绑定实际代码、完整范围和来源文件，不能只汇总成功数量。资源报告同时返回采样覆盖、实际最大间隔、SUM 的含义和中断缺口，性能优化应先定位真实分配热点，再比较冻结科学数组。

科学流程后续可以预先声明 ERCC/eGFP 技术 feature 的消融、donor 覆盖与合格的跨供体检验。当前新数据的分数一旦打开便属于已知结果，再用它调整参数只能是探索；新的确认需要另一份在评分前冻结的数据和协议。上述是后续建议，未计入本轮已实现能力或已通过验证。

## 可核查的返回入口

- [完整验收索引](F:/EACBP-architecture-p4/verification/p4-acceptance.json)：工程与科学结论分开，绑定实际代码、完整范围、参考来源、历史文件、CI XML、图与失败版本。
- [完整 native v2 独立核验](F:/EACBP-architecture-p4/verification/p4-independent-native-stability-v2.json)及[保留的 v1 失败](F:/EACBP-architecture-p4/verification/p4-independent-native-stability-v1.json)。
- [原生执行后、外部评分前的验证补充声明](F:/EACBP-architecture-p4/verification/p4-post-native-validation-addendum-v1.json)。
- [真实故障与恢复证据](F:/EACBP-architecture-p4/verification/p4-independent-recovery-v1.json)。
- [双字段外部科学评分与所有类别分母](F:/EACBP-architecture-p4/verification/p4-new-external-science-v2.json)。
- [聚类敏感性图](F:/EACBP-architecture-p4/figures/p4-clustering-stability-v1.png)与[资源观测图](F:/EACBP-architecture-p4/figures/p4-process-rss-observations-v1.png)：另有 SVG/PDF 导出，两张 PNG 已目视验收；向量版本未独立打开检查。
- [最终索引逐文件核验](F:/EACBP-architecture-p4/verification/p4-acceptance-index-verification-v1.json)：用于字节一致性，不能替代前述科学判定。

本轮能力边界为该本机 Windows 环境、冻结源码与两份数据。没有选取 winner、替换 primary 或更改科研默认；没有把完整工程验证写成生物真值已经验证。
