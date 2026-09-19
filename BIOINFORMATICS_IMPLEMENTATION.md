# 生信分析扩展实施与验收

用户授权：按架构评估提出的顺序实施，由 Luna max 负责能力实现，主代理负责集成和验收。当前工作区已有修改必须保留。

## 实施顺序

1. 必需审计项强制执行、别名明确映射；默认能力统一装配，保留注入实现，重复注册显式覆盖。
2. 真实 PyDESeq2 供体级计数模型（配对和协变量）、真实 decoupler 通路/TF 分析、供体留一法稳健性。
3. Scrublet、显式本地参考 CellTypist、满足输入条件时的 CellBender；保留原始计数和注释不确定性。
4. LIANA+ 通讯分支，明确数据库版本、样本层级和预测性质。
5. 将上述能力接入计划、路由、CLI、独立审计、证据和报告；记录同源证据依赖，显示非显著结果与敏感性边界。

## 验收要求

- 新方法必须调用真实库/可执行程序，缺依赖或输入不适用时明确失败或说明跳过原因，无模拟回退。
- 对比、原始计数、供体、配对、协变量、设计秩和参考物种等前提明确校验。
- 输出经过版本化产物和审计后才能成为证据；必需检查不能靠能力自报通过。
- 配对留一法删除完整供体，样本不足不能声称稳健；非显著和未估计值有明确状态。
- 资源文件、模型和软件版本进入可追溯记录及恢复配置；已有原始数据不变。
- 同一数据派生的分析不计为独立验证，通讯/活性为推断，不生成未经支持的因果结论。
- 针对性回归、真实库小数据运行、完整已有测试和端到端报告验证。
- 真实测序/模型数据未提供的外部程序验证边界单独记录，不以 mock 测试宣称真实生物学有效性。

## 验收状态（2026-09-18）

五步实施及软件验收已完成，编程由 Luna max 执行，主代理完成集成审查和独立验证。最终冻结回归命令：`.venv/Scripts/python.exe -m pytest -q --basetemp=temp_bio_completed_acceptance`，结果 **211 passed、414 warnings、151.12 秒**，无失败或跳过。最终公式复核发现的功能活动 OLS 置信区间问题已修复：区间与 P 值使用相同残差自由度的 Student t 分布，独立审计拒绝旧正态区间及缺失/非法自由度。完整回归包含该修复的真实库和对抗测试。

| 原始要求 | 已检查的实现及直接证据 |
| --- | --- |
| 强制审计、默认装配与注入实现 | `ScientificAuditor` 汇总独立检查并拒绝缺失必需项；`test_required_audit_coverage.py`、`test_registry_injection.py` 覆盖别名、缺失检查和显式覆盖 |
| PyDESeq2 供体计数、配对与协变量 | 真实库拟合及输入计数/设计独立重建；`test_advanced_statistics.py`，另独立探针拒绝秩亏、组内协变量变化及配对集合不一致 |
| decoupler 通路/TF 活动 | 真实 ULM、供体设计 OLS、资源版本/物种/哈希检查及 `test_advanced_pipeline.py` 报告/resume；配对/协变量真实测试验证 t 区间，篡改为正态区间或非法自由度被拒绝 |
| 完整供体留一与不确定性 | 真实三供体 refit、逐基因 summary、低样本显式跳过；完整配对流程实际完成 3 次留一拟合、总审计、汇总报告与 resume（`outputs/paired_loo_pipeline_verified/receipt.json`）；敏感性不计独立支持 |
| Scrublet 与本地 CellTypist | `test_qc_pipeline.py` 真实联合流程、未知/冲突标签、counts/ID 保留、中立证据与报告；模型不自动下载，过滤前评分篡改会被独立审计拒绝 |
| CellBender 与校正层接线 | 显式未过滤输入、可执行文件和输出；input/executable/checkpoint 固定；真实 adapter 反向子集映射审计通过（`outputs/cellbender_smoke/adapter_root_v2.json`）；完整 7 任务流程及 resume 成功（`pipeline_root_v1.json`），独立数值核对 normalization 使用 corrected_counts、原始 counts 不变 |
| LIANA 通讯与可评估范围 | 真实四供体/两条件执行、独立 BH 校验、配对与退化统计、报告/resume；缺失基因资源对及原因记录并复算；物种不匹配和无效空间坐标拒绝；真实缺失基因 fixture 通过 |
| CLI、证据与报告 | 标准 CLI 的 9 任务合成集成运行成功；同根来源追溯及保守合并；非显著/未估计/LOO 跳过保留；QC 中立结果不增加机制或因果支持；相关证据与端到端测试通过 |

**科学边界：**以上证明软件实现、输入前提、来源追溯和审计接线。合成数据/合成参考模型不能证明真实生物学有效性。真实 CellBender 公共小例子保留 ELBO 收敛及可能次优告警，`scientific_quality_certified=false`；公开单样本流程正确省略条件统计。默认供体计数模型仍用原始 counts；选择校正计数需要显式 counts_layer。参考注释冲突保留既有主标签，模型标签另列供显式选择。活动/通讯为推断，同源分析不视为独立验证。

配置入口见 `README.md`；参数及标签政策见 `ADVANCED_QC.md` 和 `ADVANCED_STATISTICS_INTERFACE.md`。运行记录及报告在 `outputs/cellbender_smoke`，测试用 Windows 与隔离 WSL 环境均有依赖快照。所有修改保存在当前工作区，未提交。

## 实施记录（历史）

以下为逐轮记录，其中“待完成”“仍 gated”等表述描述当时状态；当前状态以上方验收表为准。

- 已落盘并验证：默认能力集中装配、显式注册表保留、重复注册显式覆盖。`tests/test_registry_injection.py` 与现有编排测试合计 8 passed。
- 已落盘：证据记录全部上游根产物、同源证据维度内保守合并、敏感性不加独立支持分；报告标记启发式置信度与非独立验证。相关证据测试通过。
- 外部资源内容哈希已接入统计扩展计划；normalization 可显式选择输入层且保留 counts；数据审计识别显式元数据列；多目标请求不再静默截断。
- CLI `--advanced` / `analysis_extensions.functional_activity` 已接入PyDESeq2、decoupler和供体留一分析。双细胞检测、参考注释和LIANA已接入独立审计和计划器；背景去除分支仍待外部运行验收。
- 必需检查覆盖已接入总审计器：等价别名、缺失项拒绝、输入/输出复算、额外独立验证器接口。已修复此前仅声明的供体计数、留存率、HVG、批次混合、空间轮廓、动态基因和转移矩阵检查；保留告警级科学评估与错误级结构校验的区别。
- 新统计/QC模块与测试已落盘并注册。独立运行早版统计/QC和覆盖测试14 passed；其中统计测试真实调用PyDESeq2/decoupler，QC后端目前主要为mock边界测试，已要求补真实Scrublet/CellTypist验收。
- 已修复LOO拟合覆盖与逐基因方向稳定性的区分、未估计行状态、功能活动设计矩阵/协变量处理；后续需扩大审计与配对端到端覆盖。
- 三个Luna分支此前因额度限制中止；用户再次要求分发后，三个原有Luna max任务已恢复运行，分别负责统计、QC/CellBender及LIANA的剩余正确性和验收。主代理协调集成及隔离运行环境。
- 当前已落盘改动的完整回归：`.venv/Scripts/python.exe -m pytest -q --basetemp=temp_bio_extension_full` → **164 passed，135 warnings，96.56 秒**。其后新增的校正输入层原始计数保留测试单独运行 **1 passed**。这些结果不覆盖尚未实现的新算法，也不代表真实生物学验证。
- 接入必需检查后的全测：179 passed、1 resume测试失败（并发修改期间）；该resume单独重跑1 passed。待各分支文件稳定后重新做最终全测，此结果不视为完整验收通过。
- 新统计链端到端验证（含报告和resume）：`tests/test_advanced_pipeline.py` **2 passed**，覆盖PyDESeq2、decoupler、样本不足LOO跳过及模拟来源标记。实际拟合的小样本数据只是软件验证。
- Windows独立scrublet依赖annoy构建失败；检查已安装Scanpy源码确认其内置Scrublet实现，因此移除多余standalone scrublet依赖。已安装CellTypist 1.7.1、LIANA 1.10.0和scikit-image 0.26.0。真实Scrublet及QC边界测试9 passed。
- 稳定源码全回归：`temp_bio_stable_regression` **184 passed，195 warnings，127.01秒**（在其后CellTypist路径修复之前）。
- CellTypist实际本地模型测试发现Windows反斜杠路径触发远程模型发现，已改为绝对正斜杠路径；移除会静默丢弃用户参数的TypeError重试。使用真实SGD训练的合成模型验证load/annotate，并通过拒绝远程发现的测试替身检查本地性。合成训练仅作集成验证。CellTypist 1.7.1默认传统训练与当前sklearn的multi_class参数不兼容，SGD训练可运行；不把此软件测试视为参考模型质量验证。
- QC独立审计核对原始矩阵、细胞/基因对应、阈值和模型哈希；显式过滤要求完整过滤前评分表，核对准确保留集合并检验被删除细胞阈值。含评分篡改拒绝测试、真实Scrublet、真实CellTypist及资源测试，`temp_qc_integrity_final` **13 passed，21 warnings**。
- 外部资源固定排除output_path内容哈希，输出路径无需预先存在且不会因本次生成的输出改变恢复签名；输入文件仍要求存在及内容固定。
- LIANA已注册并开放计划分支，真实本地资源四供体两条件测试覆盖独立审计、条件比较BH校正、FDR篡改拒绝、报告和resume。`tests/test_advanced_communication.py` **3 passed，212 warnings**。排名仍非FDR，证据只描述推断，不能推导因果；Luna继续审查资源覆盖和退化统计情形。
- 最新稳定源码完整回归：`temp_bio_liana_regression` **189 passed，408 warnings，156.61秒**。这是三个Luna恢复修改之前的基线，后续变更须重新验收。

## 本轮集成验收清单（Luna max 并行实施）

| 要求 | 已有直接证据 | 待验收 |
| --- | --- | --- |
| 必需审计、注入注册表与重复覆盖 | 完整基线189项测试含required audit及registry injection | 最终稳定源码再回归 |
| PyDESeq2配对/协变量及设计前提 | 真实库测试、独立审计已有基础检查 | Luna statistics补独立设计秩/协变量检查与对抗测试 |
| decoupler通路/TF活动 | 真实ULM与统计链报告/resume测试 | 资源物种前提及设计边界深化 |
| 供体留一稳定性 | 真实refit单模块测试、样本不足跳过端到端 | 逐基因summary正式接线、完整配对LOO集成 |
| Scrublet/CellTypist | 真实算法/本地模型测试，过滤篡改审计拒绝 | 联合预处理、注释与报告端到端测试 |
| CellBender | 外部CLI适配器边界测试 | 真实CLI执行、按ID映射、来源独立审计、解除计划器gate |
| LIANA | 四供体真实库、审计、报告、resume、FDR篡改拒绝 | 可评估资源覆盖及缺失原因、退化方差、空间输入边界 |
| 证据与报告 | 同源依赖、统计非显著结果、通讯推断限定 | QC/注释观察摘要与不确定性、LOO summary展示 |

CellBender运行环境验证由主代理负责：Windows主环境不改动；WSL持久隔离路径为`/home/qiaokelihahaj/eacbp-cellbender-runtime`。CPU PyTorch已安装，CellBender 0.4.0依赖安装进行中。安装成功本身不等于算法验证，必须记录实际CLI输入、命令、退出状态、输出及审计结果。

### Luna额度中断后的检查

- 三个Luna max再次因额度限制停止；当前修改不能视为已验收。针对其落盘代码执行统计/QC/LIANA测试：18 passed、2 failed。两项失败均为LIANA执行器已新增退化方差处理，而独立审计尚未同步；已向原任务记录具体修复要求。
- WSL隔离环境已成功安装CellBender0.4.0与CPU PyTorch2.14.0，CLI帮助可运行。
- 使用官方v0.4.0 `examples/remove_background/generate_tiny_10x_dataset.py`下载公开10x mouse heart原始矩阵并裁剪，保留37760 droplets×100 genes及空液滴。目录：`/home/qiaokelihahaj/eacbp-cellbender-smoke`。
- 实际命令：`cellbender remove-background --input tiny_raw_feature_bc_matrix.h5ad --output tiny_output.h5 --expected-cells 500 --total-droplets-included 2000 --cpu-threads 4`。默认150epochs，首次运行退出码0；输入/输出形状均37760×100、cell/gene ID一致、校正值有限非负；总计数18228496→9641713。此结果验证实际运行，不代表校正生物学质量已验收。
- 首次以绝对路径调用未激活环境，附属HTML报告因找不到jupyter生成失败，尽管CLI退出码0。已用激活环境从checkpoint恢复补生成报告，结果待检查。项目adapter+审计+planner仍未完成真实联合验收，background_removal gate暂保留。
- 激活隔离环境后的checkpoint恢复已成功，退出码0并生成tiny_output_report.html。报告、指标、日志已复制到outputs/cellbender_smoke，供验收留档。
- 报告内容复核发现收敛告警：训练后半程及末epoch ELBO偏离最大值，测试末ELBO偏低，报告提示输出可能次优；指标found_cells=652（expected=500），convergence_indicator=1.233。因此真实CLI运行与矩阵格式验证通过，但不能据此声称背景校正质量通过。此告警已传给Luna QC作为审计/报告集成要求。
- 本轮Luna已修复LIANA退化统计执行/审计不一致，统计审计及LOO summary接线已落盘。冻结源码后针对统计、统计pipeline、QC、真实CellTypist、LIANA共同验收：`temp_luna_frozen_acceptance` **26 passed，255 warnings，58.60秒**。这次恢复了此前中途落盘导致的失败；不代表所有剩余能力已完成。
- WSL项目bio依赖安装完成，下一步由Luna编写独立真实CellBender adapter+auditor验证脚本。主代理同时执行完整回归。
- 最新完整冻结回归：`temp_luna_round2_full` **195 passed，390 warnings，151.31秒**。覆盖目前已保存全部源码；仍不覆盖尚未交付的QC联合流程/CellBender适配器真实联合脚本及剩余证据报告扩展，整体目标保持未完成。
- CellBender真实adapter独立复核产物：`outputs/cellbender_smoke/adapter_root_v1.json`，目标12cells×8genes反向重排，原始counts保留，ScientificAuditor结构审计通过。发现HTML/CSS `warn`误提取和AnnData警告数组序列化问题，已交新Luna小任务修复；不能将当前质量告警列表视为准确完整。
- QC/注释/background中立观察及报告指标改动已保存。现有相关回归`temp_saved_evidence_review` **20 passed，36warnings，27.41秒**；专门新增行为测试未交付，仍需补齐。
- 新批次按独立文件分工：Luna quality parser修告警解析及数组处理，Luna background plan接线并固定外部输入内容签名，Luna QC flow补真实QC/注释端到端测试。继续维持原验收范围，不因额度中断缩小目标。

### 当前补充验收（2026-09-18）

- 主代理真实执行三供体配对 PyDESeq2 留一拟合，主表与逐基因 summary 一并通过 ScientificAuditor；独立审计重建完整模型 design_rank=4、residual_df=2，并检查删除完整供体的覆盖。
- 独立设计边界探针确认：条件与 batch 完全混杂导致秩亏、同一供体条件组内协变量变化、配对供体集合不一致，均被审计输入重建拒绝。
- 标准 CLI 实际完成 Scanpy、PyDESeq2、decoupler 及留一跳过分支，9 tasks success、0 blocked。命令为 `scripts/run_standard_study.py --data outputs/advanced_cli_smoke/synthetic.h5ad --species human --tissue synthetic_test_fixture --study-id synthetic_advanced_cli --output-dir outputs/advanced_cli_smoke/runs --condition-a A --condition-b B --analysis-config outputs/advanced_cli_smoke/config.json --advanced`。报告位于 `outputs/advanced_cli_smoke/runs/synthetic_advanced_cli/61ba12bbcd2d4bbeb8fe214abb420c80/report.md`。输入是明确的合成软件测试 fixture；CLI 不从自由文本 uns 字段推断模拟来源，故该报告的数值不作生物学验证材料。
- Luna 修复了 CellBender HTML/CSS 告警误提取和 ndarray 警告展平。主代理独立运行 `tests/test_cellbender_warnings.py --basetemp=temp_warning_parser_root`：3 passed。真实保存报告中的 ELBO 偏离、wrong direction、suboptimal 与 convergence_indicator 告警得到保留。
- 仍待：CellBender 新计划合同及 checkpoint 哈希贯通真实运行；QC/注释联合流程测试；LIANA 不可评估资源对及缺失原因落入可审计产物；最后冻结全回归和文档收口。整体目标未完成。
