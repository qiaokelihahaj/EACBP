# P5：数值表合同与证据返回

本轮新增输出、临时文件、缓存、实验、日志和证据均位于 `F:/EACBP-architecture-p5`；输入沿用已有 F 盘路径并逐一绑定 SHA。代码与文档留在已有仓库。科学策略、QC 阈值和聚类参数沿用 P4；这轮验证存储与证据流程，不提出新的外部科学确认。

## 数值表的存储边界

新 TABLE 使用 `.table.npz`，一个持久文件同时包含二进制数值与 `eacbp.table.npz` / version 1 / decoder `eacbp.table.v1` 合同。合同记录矩阵形状、逐列 dtype、缺失表示、行列类型、名字与顺序。读入不启用 pickle；合同声明与数组类型、长度、键集合不符时拒绝，不通过隐式 cast 修复。

支持普通数值列、nullable 整数/浮点/布尔、Python-backed string、category、普通 Index 与 RangeIndex。原生 float32/float64 和大整数不经 CSV 类型推断。object primitive 明确区分空文本、`001` / `NA` 文本、None、NaN 和 pandas.NA。MultiIndex、datetime、complex、嵌套 object 与非空 DataFrame.attrs 目前明确拒绝；attrs 应放入产物 metadata。没有增加依赖。

旧 `.csv` 继续使用原 `pd.read_csv(path, index_col=0)` 解码；不重写历史文件、不把旧审计更新为新格式。旧 CSV 没有记录写入时的 dtype 与缺失 sentinel，因此旧源码解码相等只能证明历史读取兼容，不能补造写入类型。

差异表达结果页分页、搜索和过滤支持 typed TABLE。其展开大小在解码前检查，超出 128 MB 时目录标记不可查询；原有单文件与完整来源哈希读取预算仍适用。下载 CSV 是独立、按内容哈希缓存的文本导出，原审计文件保持原状；导出前重新核对原表哈希，缓存被改写时拒绝。完整证据包中的 typed 表保留精确类型与缺失信息。

## 正式证据验证入口

`python -m eacbp verify-evidence --index PATH --index-sha256 SHA` 只读核验被显式 SHA 绑定的 P4 acceptance envelope。返回结构化 JSON，可用 `--output NEW_PATH --report NEW_PATH` 独占发布新版本与人类可读摘要。已有输出或输入别名在验证前拒绝，不覆盖旧版本。

历史源码已演进时，另提供 `--source-archive PATH --source-archive-sha256 SHA`。归档 receipt 必须绑定同一索引与 SHA、包含唯一的原路径映射与相同的预期源码 SHA；数据证据不能借源码归档换源。原文件与显式绑定的历史副本分别报告，不宣称当前工作区仍等于旧代码。

返回分别包含工程结果、科学假设、改善条件、绝对验收、失败原因、参考字段/覆盖、来源与代码 pin、首次失败/修正补充声明/当前验证版本关系。验证读取 JSON/XML 和文件字节，不解析参考分类值，不重算 ARI/NMI，不重跑分析或授予新的科学审计。字段数量与 reported passed 不能代替完整范围、阈值、文件与版本关系核对。

退出码 0 表示证据闭合验证通过；已经记录的科学失败或指标不可估计可以得到退出码 0，其失败/null 保持原义。退出码 1 表示证据缺失、身份/范围/语义关系不成立。当前正式入口严格支持 P4 envelope；P5 工程验收索引与任意外部 JSON 不因含有 `passed` 字样获得兼容。

## 历史保护与实验方案

修改之前已逐字节归档 P4 绑定的 140 份仓库文件，包括 137 份 source map 与额外绑定文档。旧 P4 索引 SHA 为 `1282252088168836e3c9d0874e050a3960a8d2f58fd77e4c8d0bc19cd47e2662`；P5 归档 receipt SHA 为 `92dbf720384554d4d0694face74ead53c9d1179f4208f74bd0c88d3ddca5e510`。历史保护基线含 1971 个去重文件，代码副本替代已演进的当前 C 路径，旧 F 证据原位保留。

正式冻结后，分别执行：旧 10 范围/100 个输出的隔离历史加载与当前兼容读取；真实 Segerstolpe QC 105084 行表的源公式、typed 审计、跨进程与 portable 回载；PBMC 与 Segerstolpe 各一次相同科学参数的新 primary。新运行只用于核验存储变化是否影响科学输出，不读取参考标签或重算外部评分；P4 分数已知，不能称为新的预先盲确认。

P4 主科学假设仍失败：resolution 0.5 的 ARI 改善超过 0.05，但绝对 ARI 低于 0.5。次参考字段单类，分区指标不可估计。P5 不改变这些阈值、不挑选新默认参数。

## 实际验收

最终版本完整 CI 共 **1391 项通过**：core 1080、standard 240、advanced-statistics 71，0 failure/error/skip。沿用原 core 与 standard 各一项明确 deselect，没有新增排除。第一次 core 1068 项中有 6 项失败：证据描述依赖 CSV 把空文本变成 NaN，且旧测试默认新 TABLE 仍写 CSV。修正仅规范化任务/产物明确声明为空的 covariates，不把数值 NaN/Inf 或缺失估计补成有效统计量；非空声明遇空描述仍拒绝。旧 CSV 与 typed 两种格式的负例都保留。另修正 CLI 输出别名保护，覆盖不存在的声明 evidence/source/归档原文件/副本；不能把失败报告填进缺失的输入位置。

原失败 CI、日志与环境 receipt 保留。修正后依据实际第一次 CI 的包 fingerprint 精确重建并保存 133 份初版包源码，核对为 `e19752eb0fe14a25d4eba92007015ae2b3799a0c1a7235edf159234b6ee751cd`；这是事后身份匹配归档，不冒充事前冻结，也不宣称已归档初版全部测试源码。正式实验在修正后的 prospective freeze 之后启动，freeze SHA 为 `5a314f8e9886f421fcaa57b724f06444d53d51e3c3b4062f1be7c57f5a3607af`，当前包环境 source fingerprint 为 `d2c3aeb51c6452c6a321ff9f9c4b72f53a3e7cfda0ce9139d47ba744c81e8c0d`。

typed TABLE 合同/存储/审计/portable 定向回归 84 项通过；CLI/验证器初轮 62 项、修正 consumer 定向 76 项与离线行绑定 2 项通过；独立数值 helper 20 项、fresh comparator 22 项、发布 gate 54 项通过。这些与完整 CI 和各轮测试重叠，不累加为独立科研验证数量。

准备阶段通过真实轻量子进程发现 venv redirector PID 与实际 Python PID 不同；源绑定 nonce、owner 创建身份及最短父链先核验，再允许原生调用。另通过两组旧 registry 的真实 JSON 发现原稿误计 3 张 TABLE，正确为 4 张：dataset audit、QC cells、QC features 与 markers。修正保存为 bindings-v2 与独立 addendum，原稿、失败预检、初版源保持；这两处在首次原生计算前已修正。原配置除移除旧 sensitivity branch 描述与更新 study_id 外逐字段精确一致，源观察与问题合同不重新选择。额外只读收据 reader 曾把相对 metadata 路径错按 cwd 解析，失败 v1 保留，v2 按实际 artifacts 根解析；正式 native gate 自始使用正确路径。

两次真实 primary 均由当前代码完整运行，实际 launcher exit 0，各 7 个能力/10 个审计输出/4 张 typed TABLE，raw 导入另作 source pin、不制造审计。PBMC3k 保留 2562/2700 细胞、10 个簇；Segerstolpe 保留 2933/3514 细胞、30 个簇。独立科学数组与表格回归已闭合，结果见下文。

| 当前 primary | 父记录运行区间墙钟秒（含启动与结束核验） | 采样条数 | 最大观测 RSS SUM 字节 |
| :--- | ---: | ---: | ---: |
| PBMC3k | 182.492 | 217 | 780750848 |
| Segerstolpe | 377.067 | 466 | 3375902720 |

资源采样周期设为 0.25 秒，实际间隔另在 receipt 中返回。进程树观测以实际计算 Python 为根，排除启动父进程；SUM 包含重复共享页，并可能遗漏瞬时峰及未发现的后代。本轮同时有 CI 与历史字节核验的后台负载，且 P4 的资源范围包含五分支；这些数值不能构成 P4→P5 性能提升、物理峰值或生产内存限制的比较。

正式 `verify-evidence` v2 的实际退出码为 0，完整核验 1968 个去重文件身份（含索引证据与源码），明确绑定历史归档。返回 engineering passed、science failed、improvement passed、absolute acceptance failed；次参考单类、ARI/NMI null。这是原收据与文件的闭合证明，未重读参考分类值或重算外部分数。

正式独立数值兼容 v1 已通过：隔离归档旧代码与当前代码分别读取 10 个历史范围的 100 个真实审计输出，科学数组与旧表解码结果逐值、dtype、缺失与轴精确一致，容差为 0。真实 105084 行 QC 表的源公式、typed 审计、独立新进程及 portable 导入回载通过；验证前后 1971 个受保护历史文件的 SHA 一致。[数值兼容收据](F:/EACBP-architecture-p5/verification/p5-numeric-compatibility-v1.json) SHA 为 `008e1c03234f0e371ce4a034cb121af8cc444d173c08b92e675f448fc94dba0c`。这证明旧解码兼容与新类型合同保真，不能追回旧 CSV 未记录的写入 dtype 或缺失 sentinel。

该表的 mean_raw_count 与 detection_fraction 各含 52542 个有限值和 52542 个缺失值，源计数为 int64、公式商为 float64。Python float literal、round_trip 解析与源公式零容差相等；旧默认 CSV parser 分别产生 23533 个均值差异（最大绝对差 `1.8189894035458565e-12`）与 29558 个检出率差异（`1.1102230246251565e-16`），缺失位置不变。这些差异作为解析诊断返回，不转成验收容差，也不改写旧解码或历史文件。

两个新 primary 的独立回归也已完成，[fresh 比较收据](F:/EACBP-architecture-p5/verification/p5-fresh-primary-comparison-v1.json) SHA 为 `3457bcbb693782f862db5aa247f8d254adb7e57978e8e7fb3079669bcad13422`，实际退出 0，状态为 `passed_with_declared_legacy_identity_limits`。两组共 20 个审计输出中，科学数组及其 dtype、轴、保留的科学 uns 精确一致，容差为 0；仅排除事先冻结的具体来源 URI/hash/study_id 路径，没有整体丢弃 uns 或时间字段。对比器不执行分析、不读取参考标签、不重算评分，验证前后历史 1971 个文件 SHA 保持一致。

8 张表按冻结规则比较有序行列、文本值、缺失位置与数值：float64 的 Python literal 与新 typed 值精确相等；float32 则要求按新声明 dtype 解析后位精确相等，且新旧 CSV tokens 全列一致。真实数据共 508940 个有限浮点值满足 literal 精确相等，其中包含 PBMC total_counts 的 2700 个 float32 值；两类计数重叠。新表类型有明确合同，旧 CSV 原写入 dtype、轴类型与缺失 sentinel 身份仍为 unavailable，不能因表示兼容补造历史类型证明。

| TABLE | PBMC shape | Segerstolpe shape |
| :--- | ---: | ---: |
| cluster_markers | 200×17 | 600×17 |
| dataset_audit | 1×12 | 1×12 |
| qc_review_cells | 2700×9 | 3514×9 |
| qc_review_features | 130952×10 | 105084×10 |

最终[历史字节保护](F:/EACBP-architecture-p5/verification/p5-historical-protection-v1.json)已通过，逐一检查基线中 1971 个去重文件，mismatches 为空；原失败版本、P4 科学失败结果和历史归档均保留。

封存入口为 [P5 验收索引](F:/EACBP-architecture-p5/verification/p5-acceptance-v1.json)。构建器先检查冻结源码、三组 CI、数值兼容、两个真实 primary、独立回归、已闭合 CLI 与历史保护之间的实际文件、范围和版本关系，再逐文件生成证据 SHA。索引明确记录工程通过、历史科学失败以及没有新的外部科学确认；其 schema 是 `eacbp.p5.engineering_acceptance.v1`，由本轮独立构建器验证，不冒充正式 P4 CLI 支持的 envelope。

发布后，以明确的已发布索引 SHA 再逐字节检查每个新索引 evidence/source entry，结果以独立的[索引字节验证收据](F:/EACBP-architecture-p5/verification/p5-acceptance-index-verification-v1.json)为准。索引自身、其随后生成的验证收据及可能仍在写入的发布日志不加入自己的证据清单，避免循环哈希与未闭合日志；发布日志单独保留于 `F:/EACBP-architecture-p5/finalization_logs`。没有将重叠的 evidence/source entry 数量当成去重文件数。
