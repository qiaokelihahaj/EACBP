# P6：prepared 聚类敏感性实验的正式恢复入口

日期：2026-10-04–05。此阶段验证工程恢复行为，不增加科学默认参数，也不对外部生物学复现作新结论。

## 范围与存储

按用户要求，P0–P5 实验数据与计算输出已经清理，只保留代码、报告及历史验收记录。清理报告：`F:/EACBP-cleanup-2026-10-04/cleanup-report-2026-10-04.txt`；闭合核验：`F:/EACBP-cleanup-2026-10-04/cleanup-verification-v2.json`。旧验收索引反映清理前已完成的检查，当前不再具备完整产物逐字节复核条件。

P6 的新合成输入、临时文件、缓存、执行日志及证据存放于 `F:/EACBP-architecture-p6`。修改前已保存当前源码/报告快照：371 个文件、9,159,208 bytes；`verification/p6-pre-edit-snapshot-v1.json` 的 SHA-256 为 `d290aecd3383011561306fcd4c8f5c1e5ebdb4284ff6b907a6685674491d0d20`。原有工作树改动全部保留。

## 入口与信任边界

```powershell
eacbp recovery-status --prepared F:/results/prepared --bindings-sha256 <independently-checked-sha256>
eacbp resume-experiment --prepared F:/results/prepared --bindings-sha256 <independently-checked-sha256>
```

两条命令必须显式提供冻结绑定 SHA-256。状态检查返回结构化 JSON；阻断时保留原因并返回非零退出码。检查只读已有目录、绑定、冻结输入、日志及已保存的产物，不创建锁文件、不修复落后的日志 head、不补写 interrupted 事件，也不导入矩阵执行模型。合法的事件落盘而 head 尚未更新窗口可以被识别，修复留给真正的执行恢复。

诊断的 ready 只表示允许尝试恢复；不授予产物复用资格。继续执行仍通过既有原生准入流程核验任务、审计和实际产物。已有但未完整准入的研究保持原目录，需要显式进行研究级 `eacbp resume --run-dir`；不会自动删掉重算。已经有准入记录的原生目录缺失或证据改变时拒绝静默重算。源文件、冻结参数、计划、当前源码和环境必须符合原绑定。

已计划但未开始的分支可以新执行。比较与报告保存到新 invocation 版本，保留历史失败事件和此前已落盘的比较。恢复不自动选择新的聚类参数或赢家。分离的旧比较/展示文件不自动变成原生准入的可信依据；验收记录另行固定这些历史字节。

WebUI 在运行记录中提供同一诊断与恢复入口。只允许当前结果根内已有 prepared 目录，明确输入绑定 SHA，检查不提交后台工作。恢复在工作进程中再检查目录与绑定；进度读取已有 job/阶段记录，避免每次轮询重新全量哈希。未闭合阶段显示“未记录结束”，只有执行恢复可以将其追加记为 interrupted。

## 验证方法

正式实验在源码与合成输入冻结后执行：真实原生主分支完成后强制退出；第一次恢复继续剩余分支，故意在比较完成后的报告格式化阶段失败；再次恢复重新准入已有分支并生成新报告。三个执行阶段的真实子进程退出码为 23、1、0；首阶段通过 fresh API 启动，后两个阶段走正式 CLI。整个协议含 24 个子进程：1 个 fresh API、23 个 CLI，后者包括诊断和成对拒绝探测。验收交叉检查父进程记录、实际原生执行次数、准入审计、历史文件哈希及最终版本指针。

9 类拒绝实验为错误绑定、源输入改变、配置改变、计划改变、源输入缺失、主分支目录缺失、原生载荷改变、原生报告改变及并发锁。每类分别调用只读检查和执行恢复，都没有新增原生计算或新研究目录。只读诊断另检查调用前后的文件集合与哈希；单元回归覆盖落后 head 和未闭合事件。

合成数据只用于恢复与证据边界验证。工程恢复完成、任务重新准入与外部科学复现分别判断；`biological_stability_claim_supported` 仍为 false。

## 验收记录

最终受影响的 CI 分组通过 1430 项：core 1091、standard 339，失败、错误和意外跳过均为 0；各组沿用原 manifest 的 1 项明确 deselect。CI 保存实际命令、JUnit、日志、执行环境及全部代码/测试摘要，执行前后相同。未改变科学计算代码，因此没有扩大为未受影响的其他能力分组验证。恢复 API/CLI 定向回归 65 项、WebUI 定向回归 95 项均通过；这些属于 CI 覆盖的子集。F: 的独立证据 checker 另有 51 项回归通过，不能把它们与 CI 简单相加。

真实实验使用新生成的 120 cells × 40 genes 稀疏整数合成计数、固定种子和 5 个预先指定分支。v2 全程实际执行了 5 个原生研究；报告保存故障后的最终恢复重新准入并复用全部 5 个，新增原生计算为 0，已完成主分支没有重算。独立检查器核验了 35 份任务准入审计、50 个原生制品、42 个持久事件，以及 9 类成对拒绝。此前版本、失败事件和已保存比较保留。资源数值是进程身份核对后的采样 RSS 和已记录阶段耗时；强制退出阶段只有持久采样前缀，没有完整终态，不能作为物理峰值或性能基准。

WebUI 浏览器验收包含缺失目录、未开始目录和真实已完成 prepared 目录。阻断时禁用恢复按钮；通过时明确显示“尚未授权复用”，列出 5 个待重新准入分支及历史故障。检查前后 v2 的 499 个文件集合与哈希完全相同，未产生后台作业；控制台错误为 0。浏览器没有再提交恢复，以保持正式实验字节冻结；后台排队/工作进程行为由 WebUI 回归覆盖，实际恢复由正式 CLI 证明。完整页面状态为 `browser/recovery-ready-v3.ax.txt`，截图为 `browser/recovery-ready-v2.jpg`。

### 真实格式缺陷及保留的失败

第一轮 v1 在预期 ready 的检查处意外失败：产物载荷元数据实际使用 `sha256:<64hex>`，恢复 API 却沿用了仅接受裸 SHA 的文件元数据校验。原生主分支 10 个制品的字节摘要和 7 份审计实际上匹配；失败是摘要表示边界不兼容。最初 checker 的 51 项测试虽支持该真实格式，却没有调用 C: 的 inspector；C: 的 stub 使用裸摘要，因而漏过缺陷。

修复仅在 artifact payload SHA 角色接受精确、单次 `sha256:` 前缀，再执行严格小写 64 位摘要核对；绑定、输入和元数据文件的 pin 仍严格使用裸 SHA。新增 11 项真实 TABLE/JSON 序列化格式回归，覆盖改变的载荷、非法格式以及错误角色使用前缀。修复后重新冻结源码并创建独立 v2，没有改写 v1 绑定或在旧目录续跑。v1 保持失败，旧源码、原始载荷、拒绝结果和 traceback 一并保存。

经验落实为：测试 fixture 应通过真实序列化边界生成；摘要的表示与其角色必须显式区分；状态诊断与实际执行共用门禁；通过的工程验证和科学结论分别记录；修复后创建新冻结版本，保留使缺陷可复查的旧字节。

### 可复核文件与范围

以下路径均位于 `F:/EACBP-architecture-p6`，采用独立版本文件，不覆盖已有收据：

| 证据 | 相对路径 | SHA-256 |
| --- | --- | --- |
| v2 源码冻结与副本 | `verification/p6-formal-v2-source-freeze.json` | `6a9186127f5a4e5c480b1de897a9cfffa793833dc146ccd46d6b23996d8fe7c5` |
| v2 prepared 绑定 | `experiments/synthetic-recovery-v2/prepared/bindings.json` | `73d4e007a28ab326b22215016a5fa778ace8a6d6877fa63ad4eb12239dd96d84` |
| 独立真实恢复核验 | `experiments/synthetic-recovery-v2/verification/p6-independent-recovery-verification-v1.json` | `4158842726e9345f16a205f6b74522991be21d7c3ffe2b3dc736d253bdfb0d38` |
| v1 格式缺陷诊断 | `verification/p6-v1-inspector-digest-diagnostic-v1.json` | `6c968de45172eaa180a25f5fdebcebc2d756730e008aa9c307bed33aa0bd6a19` |
| 最终只读实验/CI 审计 | `verification/p6-experiment-final-audit-v1.json` | `392a7fbe5963d77e3544737743f0dec41ab52c14b1963397760e08a74019dd53` |
| 实际根进程退出观察 | `verification/p6-root-parent-observations-v1.json` | `9e67cd0bd913b55f1ccb782549e9fdbf6c915939dcf278661438ce7501fe6d3b` |
| 浏览器字节闭合 | `browser/browser-acceptance-v2.json` | `62fed2b56ce08ba8c1b0ba7b22cd3b18807d87ac377a8ae345765690e3659c46` |

实验 producer 收据保持 `verification_pending` 原字节；独立核验收据明确为 `passed`，二者不互相替代。最终原生报告位于 `experiments/synthetic-recovery-v2/prepared/comparisons/9e279d36d2874b84af5ab72e3c9eb69a/report.md`；其 SHA-256 为 `6d5a856b7d35f557e88907dba85dc425fd1ff68a3d0f145a3d2dddd65ed9ec76`。根实验驱动的退出码通过实际工具返回值记录为观察收据；实验子进程的 23/1/0 来自已闭合 parent JSON，可另行交叉核对。

最终索引为 `verification/p6-acceptance-v2.json`，独立逐字节复核结果为 `verification/p6-acceptance-index-verification-v2.json`。索引固定 F: 的源码副本、真实输入/产物、版本历史、日志、收据和报告；它证明这些保存字节及当时受影响 CI 的结果。它不独立重算模型、确认生物学结论或授予未来的 native 复用；以后改动当前源码/环境，实际恢复仍可能按严格身份规则拒绝。此前 v1 索引及其通过记录保留为补清理 WebUI 临时 fixture 之前的历史记录，其被清理文件不再可复核，当前完整字节核验使用 v2。

```powershell
& C:/Users/hahaj/Desktop/EACBP/.venv/Scripts/python.exe -B F:/EACBP-architecture-p6/verify_p6_index.py --index F:/EACBP-architecture-p6/verification/p6-acceptance-v2.json --index-sha256 <独立核对的完整索引SHA> --output F:/EACBP-architecture-p6/verification/new-index-check.json
```

索引 SHA 记录在上述独立核验结果中；复核时必须显式传入，不自动把待验文件的当前哈希当作可信来源。输出路径必须尚不存在。

本轮临时清理记录为 `verification/p6-temporary-cleanup-verification-v1.json` 与 `verification/p6-webui-temporary-cleanup-verification-v2.json`。仅清理 P6 的 `tmp`、`cache`、`test-tmp`，共 70,214 个临时文件、971,118,907 logical bytes；临时根已不存在。清理前后 1,951 个原生实验/源码副本文件、56,527,575 bytes 完全一致；两轮小型合成实验和源码/报告证据保留。新阶段没有重新下载已清理的公开数据。工程验收状态为通过；科学状态仍为 `not_evaluated`，`biological_stability_claim_supported=false`。
