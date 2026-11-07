# 公开数据验证记录

记录日期：2026-10-01。本文记录 PBMC3k 和 Kang 2018 公开数据的固定输入、当前工程运行结果和可复跑命令。

**状态：本轮工程和方法验收通过，领域审查仍待完成。** 已下载并核验公开输入，完成 PBMC3k 独立双跑、恢复和研究包往返，以及 Kang 2018 配对供体统计与独立科学审计。Windows 和 Docker/Linux 全量 CI 均通过未预期 skip 门禁；真实研究结论仍为 `pending_domain_review`。

## 固定输入与来源

输入存放在 `F:\EACBP-public-validation`。下载器 receipt 位于 `F:\EACBP-public-validation\receipt.json`，保留直连下载 URL、来源、引用、许可说明、预期与实测哈希、文件大小及核验时间。最终证据索引为 `references/acceptance-final.json`，记录 11 项关键证据的文件哈希和公开输入的最终核验哈希。

| 数据 | 固定路径 | 形状（细胞 × 特征） | SHA-256 | 来源与许可 |
|---|---|---:|---|---|
| PBMC3k raw | `datasets/pbmc3k_raw.h5ad` | 2,700 × 32,738 | `89a96f1beaa2dd83a687666d3f19a4513ac27a2a2d12581fcd77afed7ea653a1` | [Scanpy raw 数据说明](https://scanpy.readthedocs.io/en/latest/api/generated/scanpy.datasets.pbmc3k.html)；上游 [10x PBMC3k 页面](https://www.10xgenomics.com/datasets/3-k-pbm-cs-from-a-healthy-donor-1-standard-1-1-0) 明确标为 CC BY 4.0。Scanpy 序列化文件没有单独的数据文件许可声明。 |
| PBMC3k processed reference | `references/pbmc3k_processed.h5ad` | 2,638 × 1,838 | `0db367b991dd95809732b218539ede489bea99113807f62ebd7ccc970025fe38` | [Scanpy processed 数据说明](https://scanpy.readthedocs.io/en/1.12.x/generated/scanpy.datasets.pbmc3k_processed.html)；其上游 10x 来源许可为 CC BY 4.0。Scanpy 没有为这个 processed H5AD 单独声明许可；receipt 只记录上游许可，不把 Scanpy 软件许可当作数据许可。 |
| Kang 2018 | `datasets/kang_2018.h5ad` | 24,673 × 15,706 | `e6a5adac64dcdeb36eaba27db49b63e0c64bb0ed4a64c6705971506b41c39830` | [Figshare article 19397624](https://figshare.com/articles/dataset/Kang_HM_Subramaniam_M_Targ_S_Nguyen_M_et_al_2017/19397624)，version 2，file id `34464122`，CC BY 4.0；文件 MD5 `adb2246232e8493031c576982c0c02a3`。引用 Kang et al., “Multiplexed droplet single-cell RNA-sequencing using natural genetic variation,” DOI [10.1038/nbt.4042](https://doi.org/10.1038/nbt.4042)，GEO `GSE96583`。 |

PBMC raw 和 processed 文件来自同一 10x 健康供者示例；processed reference 仅用于聚类完成后的比较，没有进入数据准备、训练或聚类。PBMC3k 没有 donor/condition 设计，因此不执行差异表达或 donor-level 推断。

## 环境与隔离方式

最终 `eacbp:py312-cpu` Linux/amd64 镜像使用 Python 3.12.14，镜像 ID 为 `sha256:3cfdb1134215ee4a6b0a524ff61519f95b7f8aa5238dd81ca52ab487af58aa16`，大小 395,395,394 字节，默认 UID 为 `10001:10001`。97 个科学依赖仍由版本和下载哈希固定的 lock 安装。`references/docker-image-seeded-final.json` 与 `references/source-identity-seeded-final.json` 记录完整镜像、版本和代码指纹；镜像已安装包的 111 个 Python 文件与该阶段验收时工作树逐一字节一致。revision 标签含 `public-seeded-validation-dirty`，表明这次验收包含未提交源码。

最终真实分析只挂载 `scripts` 到 `/harness/scripts`，使用镜像内 `/opt/venv/lib/python3.12/site-packages/eacbp`，没有用工作树覆盖已安装包。原始输入只读，容器使用 `--network none` 和只读 root filesystem；真实结果及 `/scratch` 缓存、临时目录绑定到 F:。早期诊断使用过只读工作树和 `PYTHONPATH=/workspace`，其日志保留，最终验收以已安装包运行和源文件一致性检查为准。

一次镜像重建遇到 PyPI 文件读取超时。Dockerfile 将 pip 下载超时设为 120 秒，把锁定科学依赖安装放在应用 wheel COPY 之前，并将 revision ARG/LABEL 移到末尾。完整重建成功；随后只改变 revision 的构建中科学依赖和应用安装层均命中缓存。证据为 `runs/docker-cache-final-build.log`、`runs/docker-revision-cache-probe.log`。最终镜像的默认 scratch/non-root 探针、只读运行探针及合成 CLI 生命周期均通过，见 `runs/container-seeded-final.log`；合成输出仍明确标记 `is_simulated=true`。

## PBMC3k 最终独立双跑

最终两个独立容器的输出目录为：

- `runs/pbmc-seeded-a/pbmc3k-validation-acd845754ce244039d437aa675393f82`
- `runs/pbmc-seeded-b/pbmc3k-validation-b74e372ca2ac40488397bd79776e6fbf`

两次均为 5 个成功任务、0 个阻塞任务；QC 从 2,700 个细胞保留 2,699 个，保留 32,738 个特征，得到 11 个 cluster。源文件前后 SHA 相同。`runs/pbmc-seeded-reproducibility.json` 确认细胞/基因顺序、counts、PCA、UMAP、邻接距离与连通图、聚类和 cell-type 标签全部相同；PCA、UMAP 和二维 embedding 按位一致。

在相同最终代码下，run A 的 `resume`、`report`、`export`、`import`、导入后 `report` 全部通过，导入前后 Markdown 报告相同，原始数据未改变。证据为 `runs/pbmc-seeded-lifecycle-final/lifecycle-result.json`；研究包大小 69,242,043 字节，SHA-256 为 `f392313592baf1fddbf97c663df31159bb7b2b440b29fab2efd679be0b5272a5`。

在查看参考标签前固定了标准分析配置：`min_genes=10`、`max_mito_pct=20`、`no_correction_v1` integration、`scanpy_leiden_umap_v1` clustering，以及下列 marker 集。QC 从原始计数计算线粒体比例；本次识别到 13 个线粒体基因。固定 marker 用于分析后的粗粒度 cell-type 注释，不由参考标签调整。

| 注释 | 预先固定的 marker |
|---|---|
| B cells | `MS4A1`, `CD79A`, `CD79B`, `CD74` |
| CD4 T cells | `IL7R`, `LTB`, `CCR7`, `CD3D`, `CD3E` |
| CD8 T cells | `CD8A`, `CD8B`, `CCL5`, `GZMK`, `CD3D` |
| NK cells | `NKG7`, `GNLY`, `PRF1`, `KLRD1`, `GZMB` |
| CD14 monocytes | `LYZ`, `S100A8`, `S100A9`, `FCN1`, `CTSS` |
| FCGR3A monocytes | `FCGR3A`, `MS4A7`, `LST1`, `IFITM3` |
| Dendritic cells | `FCER1A`, `CD1C`, `CLEC10A`, `CST3` |
| Platelets | `PPBP`, `PF4`, `SDPR` |

聚类完成后，脚本才打开 processed reference，并按精确 barcode 交集对齐：参考有 2,638 个细胞，聚类结果有 2,699 个细胞，交集为 2,638 个；多出的 61 个聚类细胞没有参考对应项。对 `louvain` 的聚类比较得到 ARI `0.5376852037`、NMI `0.7260149863`。按固定的人 PBMC 标签别名评估 marker 注释，2,638 个共享细胞上的 accuracy 为 `0.9298711145`、macro-F1 为 `0.8686562326`、ARI 为 `0.8623993998`。

这些是这次数据和固定流程的描述性指标，不是“所有细胞类型都可靠”的判定。特别是参考中的 37 个 Dendritic cells 只有 20 个与预测标签一致；稀有类的可信度和生物学意义仍需领域审查。run A 保存 `public_pbmc_validation.json`、`annotation_evaluation.json` 和 `annotation_confusion.png`。

历史目录 `runs/pbmc3k-validation-cc0c35d8fbe54fbe8e0d9789da4449c5` 与 `runs/pbmc-installed-final` 均保留。旧 run 在后续源码修改后 resume 正确拒绝不同的 source fingerprint。两次早期分析还暴露了 PCA 上游随机性：输入、QC 和 normalized 制品相同，但 PCA 和图发生变化，11/10 个簇导致注释指标不同。`references/pbmc-prior-reproducibility-probe.json` 保存数值比较。修复采用 contract seed 生成稀疏 ARPACK 初始向量；没有根据参考标签调参。SciPy 的初始向量定义见 [svds 文档](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.svds.html)。

## Kang 2018 最终工程复跑

该研究包含 8 位 donor 的 `ctrl`/`stim` 配对样本。中间运行对 CD14+ Monocytes（5,697 个细胞）及 CD4 T cells（11,238 个细胞）分别保留全部 15,706 个基因，形成每个对比 16 个 donor-condition pseudobulk 样本；模型为 `~condition + donor`，对比为 `stim` versus `ctrl`。donor 是统计单位，计算结果不能替代对样本设计和 donor 独立性的判断。

最终已安装镜像的真实复跑位于 `runs/kang-cd14-seeded-final` 和 `runs/kang-cd4-seeded-final`。summary 分别报告 `overall_status=passed`，engineering、method acceptance、independent scientific audit 及七个预先指定基因的方向检查全部通过，输入 source hash 未变化。CD14+ Monocytes 子集为 5,697 × 15,706，CD4 T cells 子集为 11,238 × 15,706；每项比较均使用 8 位 donor 的 16 个配对 donor-condition pseudobulk 样本和 `~condition + donor` 模型。CD14 与直接 PyDESeq2 参考调用的效应/区间最大差异为 `1.776e-15`、FDR 最大差异为 `1.110e-16`，包括不可估计项的掩码一致性检查。

固定方向检查基因为 `ISG15`, `IFI6`, `LY6E`, `MX1`, `IFIT1`, `IFIT2`, `ISG20`；来源为 [Seurat Kang 教程](https://satijalab.org/seurat/articles/integration_introduction.html)。这里只比较预先指定的正向效应，不要求 donor-level p/q 值与教程的 cell-level 测试一致。`runs/kang-marker-effects.png` 和 `.pdf` 显示效应及 95% Wald 区间。此前发现的 `label` 条件别名现映射到建模 `condition`，独立审计重建 contrast，并仍拒绝真实篡改。

这些结果通过工程和方法审计，但仍不能代替领域判断。当前研究验收状态为 `pending_domain_review`，`biological_validity_claimed=false`。

## 缺陷修复与回归状态

公开数据流程暴露并促成了几项修复：QC 根据原始 counts 计算线粒体比例，无法识别线粒体基因时不把缺失的 `percent_mito` 当作 0；人类 PBMC marker reference 不再被脑组织默认值覆盖，空 marker 集不能给出假阳性；acceptance runner 的 `_filled` 名称错误已修正；contrast 执行和独立审计都解析实际模型 factor，并保留篡改拒绝检查；Integration/Subset 的 contract seed 现在确实传入稀疏 PCA。公开脚本验证原始文件不变、拒绝覆盖已有输出，并用跨平台纯路径测试验证绝对路径、外部目录和 C: 拒绝策略，未放宽 CI skip 门禁。

| 最终检查 | 结果 | F: 证据（相对验证根目录） |
|---|---|---|
| Windows Python 3.13 全量 CI integration | 587 passed，1 expected skip；384.34 秒；JUnit 门禁通过 | `runs/public-ci-seeded-final.log`、`.xml` |
| Docker/Linux Python 3.12.14 全量 CI integration | 586 passed，2 expected skips；175.14 秒；JUnit 门禁通过 | `runs/public-ci-linux-ram-final.log`、`.xml` |
| 实际稀疏 PCA、Integration/Subset 回归 | 10 passed | `runs/pca-regression.log`、`.xml` |
| 受影响能力、公开脚本和真实方法定向回归 | 29 passed，0 skipped | `runs/public-pca-final.log`、`.xml` |
| 跨平台公开脚本与 CI 分组定向检查 | 27 passed，0 skipped | `runs/public-path-final.log`、`.xml` |

两份全量报告都包含 588 项测试，没有未预期 skip。共同跳过项是未提供外部 CellBender 真报告；Linux 另跳过 Windows-only job-assignment 测试。定向检查覆盖部分相同测试，不与全量计数相加。

Linux 测试在最终相同代码、97-pin 科学环境中执行。额外三项测试工具固定为 pytest 8.4.2、iniconfig 2.3.0、pluggy 1.6.0，安装到 F:，没有覆盖科学依赖，安装 receipt 为 `references/docker-test-deps-installed.json`。Windows 共享目录上的第一次 FD 捕获失败、随后共享目录复跑很慢，保留日志后停止；最终将经逐文件核验的代码快照和合成测试临时文件放在 RAM tmpfs，公开数据、库缓存、日志和 JUnit 仍在 F:。驱动脚本为 `references/run_linux_ram_ci.py`。

首轮历史全量快照 `runs/public-full-regression.log` 为 578 passed、2 failed、1 skipped。失败分别来自新测试未分配 CI lane，以及源码在测试期间变化触发正确的严格恢复拒绝；后者无需放宽恢复签名。后续六项审计回归失败来自旧测试夹具缺少实际 formula，已补上 `~condition` 并保留全部篡改断言。所有中间失败和最终结果均保留。

## 可复跑命令

在仓库根目录用 PowerShell 执行。每次复跑都建立新的输出目录：PBMC 脚本会在其下创建唯一 run 子目录；Kang 脚本拒绝覆盖已有输出。数据和源码只读，F: 的 `cache` 目录挂到 `/scratch`，临时文件随镜像的 `TMPDIR=/scratch/tmp` 落在 F: 缓存目录中。镜像的非 root 用户为 `10001:10001`。

```powershell
$repo = (Resolve-Path .).Path
$repo = $repo.Replace('\', '/')
$dataRoot = 'F:/EACBP-public-validation'
$scratch = "$dataRoot/cache"
$tag = (Get-Date -Format 'yyyyMMdd-HHmmss') + '-' + [guid]::NewGuid().ToString('N')
$pbmcWork = "$dataRoot/runs/replay-$tag/pbmc"
$kangCd14Work = "$dataRoot/runs/replay-$tag/kang-cd14"
$kangCd4Work = "$dataRoot/runs/replay-$tag/kang-cd4"
New-Item -ItemType Directory -Path $pbmcWork, $kangCd14Work, $kangCd4Work | Out-Null

$dockerArgs = @(
  '--rm', '--network', 'none', '--read-only', '--user', '10001:10001',
  '--mount', "type=bind,source=$repo/scripts,target=/harness/scripts,readonly",
  '--mount', "type=bind,source=$dataRoot,target=/validation,readonly",
  '--mount', "type=bind,source=$scratch,target=/scratch",
  '-e', 'HOME=/scratch/home',
  '-e', 'XDG_CACHE_HOME=/scratch/cache',
  '-e', 'NUMBA_CACHE_DIR=/scratch/cache/numba',
  '-e', 'MPLCONFIGDIR=/scratch/cache/matplotlib',
  '-e', 'TMPDIR=/scratch/tmp'
)

docker run @dockerArgs --mount "type=bind,source=$pbmcWork,target=/outputs" `
  eacbp:py312-cpu python /harness/scripts/verify_public_pbmc.py `
  --data /validation/datasets/pbmc3k_raw.h5ad `
  --reference /validation/references/pbmc3k_processed.h5ad --work-dir /outputs

docker run @dockerArgs --mount "type=bind,source=$kangCd14Work,target=/outputs" `
  eacbp:py312-cpu python /harness/scripts/verify_public_kang.py `
  --data /validation/datasets/kang_2018.h5ad --work-dir /outputs `
  --cell-type 'CD14+ Monocytes' `
  --lockfile /opt/eacbp-build/requirements-scientific.lock

docker run @dockerArgs --mount "type=bind,source=$kangCd4Work,target=/outputs" `
  eacbp:py312-cpu python /harness/scripts/verify_public_kang.py `
  --data /validation/datasets/kang_2018.h5ad --work-dir /outputs `
  --cell-type 'CD4 T cells' `
  --lockfile /opt/eacbp-build/requirements-scientific.lock
```

全量 CI 回归需要已安装完整开发与科学测试依赖的项目 `.venv`；此运行器执行 integration suite，覆盖所有测试并用 JUnit 门禁检查未预期 skip。临时目录和报告都指定到 F:：

```powershell
$testWork = "F:/EACBP-public-validation/runs/replay-$tag/tests"
New-Item -ItemType Directory -Path $testWork, "$testWork/cache/numba", "$testWork/cache/matplotlib", "$testWork/tmp" | Out-Null
$env:TEMP = $testWork
$env:TMP = $testWork
$env:TMPDIR = $testWork
$env:XDG_CACHE_HOME = "$testWork/cache"
$env:NUMBA_CACHE_DIR = "$testWork/cache/numba"
$env:MPLCONFIGDIR = "$testWork/cache/matplotlib"
$env:PYTHONDONTWRITEBYTECODE = '1'
$env:PYTEST_ADDOPTS = '-p no:cacheprovider'
.\.venv\Scripts\python.exe -B scripts/verify_ci_test_suite.py integration `
  --junitxml="$testWork/public-full-regression.xml"
```

本轮仅覆盖上述公开 PBMC/IFN 分析及工程恢复；PBMC marker 注释和 Kang 研究解释仍需领域审查。Slurm/Apptainer 目标集群、真实 CellBender 报告以及其他组织、疾病和方法的研究适用性不由本轮验收确认。

## 2026-10-01 追加复核与架构经验

继续使用 F: 上已有的固定公开输入，没有重复下载或修改原始数据。追加验收记录为 [final-review.json](F:/EACBP-public-validation/runs/goal-final-review-20261001/final-review.json)，保留新运行日志、JUnit、源码哈希和数值比较；上面的历史验收文件没有被覆盖。三份输入和两份既有验收索引中的 **38 项哈希引用**全部实测匹配，引用数包含重复引用，不代表 38 个独立数据集。

本次发现并修复了 PBMC 基准入口的验收缺口：原脚本记录任意 raw 文件哈希，却没有强制要求它与下载目录固定文件相同；processed reference 没有固定哈希；流程成功即可返回总状态 success，即使参考比较只覆盖部分细胞。当前 `verify_public_pbmc.py` 从下载 catalog 读取两份固定 SHA-256，运行前核验、运行结束再次核验两份输入，并分别保存 `statuses.engineering` 与 `statuses.benchmark`。只有工程成功、全部参考 barcode 有有效聚类标签、全部参考标签能按固定别名评分且两份输入未改变时，总状态才为 success。输入替换、运行中参考变更或覆盖不完整会返回失败 exit code。ARI、accuracy、F1 没有加入依据已见结果设定的通过阈值；benchmark passed 表示输入与比较范围核验完成，科学解释仍为 `pending_domain_review`。

追加定向测试为 **39 passed、0 skipped**，覆盖固定输入、参考替换、运行中变更、部分匹配、缺失/未知参考标签、原始 counts QC、marker 审计和 CI 分组。日志和 JUnit 位于 `runs/goal-final-review-20261001/focused-tests.log` 与 `.xml`。这些是本次定向测试，不与上面的全量计数相加。

随后用 `eacbp:research-value-cpu` 已安装包进行新 PBMC 真实复跑；镜像的 111 个 EACBP Python 文件已重新核验与当前源码一致。数据和 harness 只读挂载，禁用容器网络，新实验文件、缓存和临时目录均落在 F:。新运行位于 `runs/goal-final-review-20261001/pbmc/pbmc3k-validation-68afe084757c434987fb074f15d0b619`；engineering 和 benchmark 均 passed，2,638/2,638 个参考细胞全部完成比较。新 H5AD 与先前 seeded run A 的 11 项检查完全一致，包括 counts、细胞/基因顺序、PCA、UMAP、二维 embedding、两个邻接矩阵以及 cluster、Leiden、cell-type 标签。注释 accuracy 仍为 2,453/2,638（92.99%），macro-F1 为 0.868656；聚类 ARI 为 0.537685、NMI 为 0.726015。

独立只读复核还从既有 Kang 主表重新计算了 marker 与 donor 敏感性结果，并对保存的直接参考表再次逐基因比较：两个细胞类型各 15,706 个基因、9 个拟合、5 项统计量和 NA/Inf mask 均一致（`rtol=1e-6`、`atol=1e-8`，最大绝对差约 `3.55e-15`）。没有重新拟合 reference；具体候选与机器性能范围仍以 [RESEARCH_VALUE_VALIDATION.md](RESEARCH_VALUE_VALIDATION.md) 为准。

这轮复现沉淀的架构经验如下：

| 实验暴露的问题 | 已落实的处理 | 适用边界 |
|---|---|---|
| 文件能读入并不保证是固定基准输入 | 下载和执行共享 checksum catalog；两输入前后核验 | 本地哈希相对固定记录验证身份，不是外部数字签名 |
| 工程成功可能掩盖比较范围不完整 | 工程状态、完整参考比较、科学解释分别记录 | 完整比较不等于准确率或生物学有效性通过 |
| 稀疏 PCA 的随机初始化影响双跑结果 | contract seed 传入稀疏求解器，并实际比较 PCA、图和标签 | 当前固定输入/环境下的可重复性，不保证跨库版本按位一致 |
| 表达变换可能污染 QC 和统计输入 | QC 使用原始 counts；pseudobulk 保留全基因和合法 NA | donor/condition 元数据仍需描述真实实验设计 |
| 供体剔除后旧模型计数会误导审查 | 每个 fold 重建样本数、donor 数、design rank 与独立审计 | LOO 是选择后敏感性检查，不能替代独立队列 |
| 总体注释指标掩盖稀有细胞表现 | 保存逐类混淆表，并明确 DC 只有 20/37 正确 | 不从约 93% 总体准确率推断全部细胞类型可靠 |

科学结论限于固定 PBMC 参考比较和 Kang 预设 IFN 响应方向的复现。PBMC 的 QC/Leiden 流程不同于历史 processed reference 的原流程；Kang 原论文还包含 demuxlet 和遗传分析，本轮没有复现整篇论文。直接 PyDESeq2 数值比较及 LOO 都使用同一数据，不能据此宣称独立生物学验证或新发现。
