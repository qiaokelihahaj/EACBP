# EACBP 容器化规划

日期：2026-09-19。状态：按当前代码重新核实并更新规划，尚未实现或构建镜像。代码核查基线：`a11f905`；本次规划文档修改不属于该提交。

## 1. 决策摘要

采用 **Linux amd64、CPU 批处理镜像**。本地及 CI 使用 Docker；现有 Slurm 集群优先评估 Apptainer，以同一 OCI 镜像转换为 SIF 执行。一次容器执行对应一个研究作业，结果落到持久化目录。

先交付 h5ad 分析镜像及同镜像恢复验收，再增加高级分析、FASTQ 工具和集群支持，最后完成发布与跨节点恢复评估。现有代码没有 Web 服务、数据库或消息队列，首期无需端口、Compose 服务编排或 Kubernetes。标准分析采用 CPU；可选 CellBender 的 CPU/GPU 运行时单独评估。Slurm 使用名为 gpu 的分区不等于已经申请 GPU。

本次修订纠正旧规划的两个前提：项目已有可安装的 `eacbp` CLI，也已有 `resume` 命令。容器化应复用这些入口，而不是重新开发一套运行与恢复接口。

规划默认保留现有计算方法、审计、证据和失败语义；容器化的成功标准是环境与执行可复现，不是改变科学结论。

## 2. 已核实的项目现状

| 代码或文件 | 观察 | 对容器化的影响 |
|---|---|---|
| `pyproject.toml` | Python >=3.10；另含 `advanced-statistics`、`advanced-qc`、`communication` 等扩展；多数依赖仅设下限 | 必须按镜像功能独立生成 Linux 依赖锁，不能用一次普通 pip 安装作为发布依据 |
| `environment.yml` | Python <3.12，未完整覆盖 standard/fate | 属于旧环境入口，应明确适用范围，不能作为新镜像唯一依赖来源 |
| `requirements-tested-windows-py313.txt` | Windows/Python 3.13 快照，包含 editable Git SSH 安装 | 只作版本参考；禁止直接拿来构建 Linux 发布镜像 |
| `pyproject.toml`、`eacbp/cli.py` | 已声明 `eacbp = eacbp.cli:main`；含 plan/run/resume/inspect/cleanup/report | 安装 wheel 即可提供主入口，无需复制脚本来支持 h5ad |
| `eacbp/cli.py` | 新建使用 `--run-dir`，拒绝已存在目录；恢复读取保存的配置与已导入产物 | 容器首期即可验收运行与恢复；不要提前创建具体 run 目录 |
| `scripts/run_standard_study.py` | 旧 h5ad 脚本仍保留 | 可作为兼容入口，不承担新容器的恢复契约 |
| `scripts/run_fastq_to_biology_pipeline.py` | real/demo 显式区分；调用 kb 或 STAR；返回非零失败状态 | FASTQ 镜像必须包含实际工具链，不能用 demo 证明真实比对可用 |
| `eacbp/artifact/` | 使用文件锁、硬链接、原子替换与事务目录 | 持久化卷必须验证这些文件系统语义，不能直接替换成对象存储挂载 |
| `eacbp/orchestrator/checkpoint.py` | 恢复指纹含源代码、依赖和 `platform.platform()` | 镜像相同也不保证跨宿主内核恢复；需要专项验收 |
| `slurm/*.sbatch` | 宿主 Python、项目和日志路径存在绝对路径默认值 | 新增容器模板，保留现有运行方式作为迁移回退 |
| `.github/workflows/tests.yml` | 已配置 Windows/Linux、Python 3.11/3.12/3.13 测试和 wheel 构建 | 可扩展容器测试；配置存在不代表远端测试已成功 |

本次只读探测发现 `docker.exe` 与 `wsl.exe`，但 `docker version` 无法连接 `dockerDesktopLinuxEngine` 命名管道，当前会话尚不具备可用的本地 Docker Linux 构建服务；未启动或修改 Docker/WSL。未连接集群核查 Apptainer、节点权限、存储与配额。核查开始时已跟踪文件无修改，但存在大量未跟踪的 `temp_*` 验证目录，必须从构建上下文排除。实施时仍需选定可追溯源码快照。

## 3. 镜像边界与依赖策略

### 3.1 发布目标

| 目标 | 内容 | 交付阶段 |
|---|---|---|
| `analysis` | EACBP + standard + fate；支持 h5ad、空间分析和现有本地知识模块 | P0 |
| `analysis-test` | 与 analysis 同一运行依赖，加测试工具和测试文件 | P0，仅用于验证 |
| `analysis-extended` | analysis + advanced-statistics + advanced-qc + communication；本地模型与网络资源外部挂载 | P1，高级功能有独立验收，不把依赖缺失的跳过算通过 |
| `fastq-kb` | analysis + 固定版本 kb-python/kallisto/bustools | P1，优先对应当前默认路径 |
| `fastq-star` | analysis + 固定版本 STAR 及 gzip 解压所需命令 | P1，按 STARsolo 用户需求验收 |
| `cellbender` | analysis-extended + 独立 `/opt/cellbender` 环境，执行适配器指定的绝对路径 | 可选 P1；CPU 先验收，GPU 版本按实际需求另行锁定与验证 |

不先拆分编排器、审计器和能力模块为多个服务：它们当前共享本地事务注册表与进程内状态，拆分会增加新的协议和一致性工作。

### 3.2 构建约定

1. 候选基线为 Python 3.11 + Debian bookworm slim，目标 `linux/amd64`；这是兼顾现有环境约束的待验证选择。若当前依赖无法解析，应评估调整版本或 Python 3.12，并同步迁移文档，不能静默降级科学方法。
2. 在目标 Linux/Python 环境解析完整传递依赖并保存精确版本和哈希；保留 `harmonypy==0.0.10`、CellRank/AnnData 范围等现有约束。运行与测试使用相同科学计算依赖锁；高级扩展联合解析，不在基础镜像上无约束追加 pip 安装。构建工具也固定版本，离线安装关闭隐式构建依赖下载。
3. 基础镜像固定实际可获取的 digest；Python 包构建 wheelhouse，应用安装非 editable wheel。基础系统包与外部工具记录版本、来源、校验和；不使用浮动最新版工具下载。
4. 多阶段构建：builder 负责 wheel/必要编译；runtime 仅保留运行依赖。所需动态库通过安装与导入测试确定，不预设一长串未验证系统包。
5. 主入口采用安装 wheel 提供的 `/opt/venv/bin/eacbp`，支持 `python -m eacbp`。只有仍使用旧 FASTQ/兼容入口的 target 才显式复制所需脚本至 `/opt/eacbp/scripts/`。包与虚拟环境放在 `/opt` 下，不依赖用户 home。
6. 禁止 `COPY . .` 无差别打包。应用源码、必要元数据和入口采用白名单 COPY；`.dockerignore` 排除 `.git`、`.agents`、`.venv`、`.artifacts`、所有 pytest/temp 目录、outputs/logs/build、原始数据与本地配置。
7. 镜像记录源码 revision、源码内容哈希、依赖锁哈希、基础镜像 digest 与构建时间。发布镜像从确定的源码快照构建；开发脏工作区构建必须有独立标识。
8. kb 安装后验证实际 `kb`、`kallisto`、`bustools` 可执行且运行时无需下载；STAR 验证版本、动态链接和压缩输入路径。工具版本要与引用索引制作版本一起归档。
9. CellTypist 模型、decoupler 网络、LIANA 的用户资源、参考索引和 CellBender checkpoint 在部署前准备并记录版本/哈希。选定离线资源的验收作业关闭网络运行，以发现隐式下载。模型分发权限按所选资源核实，不把用户数据或参考模型默认烘焙进镜像。
10. CellBender 适配器在当前容器内直接调用外部进程，因此其可执行文件及依赖必须在该容器可见；不能仅启动另一个容器就认为适配器能够访问它。独立 Python 环境隔离其依赖，仍在同一作业镜像执行。GPU 变体需另验 NVIDIA 驱动、容器运行时、PyTorch/CUDA 组合与 Slurm GPU 申请，基础镜像不承担该保证。

Docker 官方建议使用多阶段构建、固定基础镜像 digest 和非 root 用户，参见 [Docker 构建最佳实践](https://docs.docker.com/build/building/best-practices/)。具体科学计算依赖的兼容性仍以本项目 Linux 实测为准。

## 4. 文件系统与执行契约

| 容器路径 | 权限 | 内容与生命周期 |
|---|---|---|
| `/opt/eacbp`、`/opt/venv` | 只读 | 脚本、应用与依赖 |
| `/data` | 只读绑定 | FASTQ/h5ad 原始输入，manifest 内使用容器可见路径 |
| `/refs` | 只读绑定 | kb 索引、t2g、STAR genomeDir、whitelist、GTF、本地模型和分析网络 |
| `/config` | 只读绑定 | 样本 manifest、markers、terminal states 等运行配置 |
| `/outputs` | 持久化读写绑定 | 完整 run 目录、注册表、事务载荷、journal、报告和比对结果 |
| `/scratch` | 作业级读写绑定 | 可丢弃的缓存和临时文件；使用节点本地磁盘时不提供恢复保证 |

固定容器内路径，允许宿主目录变化。样本 manifest 中的宿主绝对路径不能原样传入；准备容器版本 manifest 并保留原文件。所有输入路径在运行前检查可读，输出路径检查可写。

**必须保存整个 run 目录。** `_transactions` 中包含已经提交的真实 artifact，不可视为普通缓存批量删除。索引、载荷与 journal 应在同一个支持所需原子操作的文件系统中；快照或备份在作业停止/一致性点进行。

存储上线前在实际挂载点运行锁、硬链接、原子替换和多进程冲突测试。NFS/并行文件系统、Windows 绑定目录均不能只凭名称假设兼容；不满足时先使用经验证的 Linux 存储。每个研究使用独立 run 目录，不把当前注册表当作分布式数据库。

Docker 使用非 root 用户并支持宿主 UID/GID 对齐；输出目录由部署者预先准备。默认工作目录 `/outputs`，所有脚本仍显式接收 `--output-dir`。建议 `PYTHONUNBUFFERED=1`、`PYTHONDONTWRITEBYTECODE=1`、`PYTHONNOUSERSITE=1`，并把 `HOME`、`XDG_CACHE_HOME`、`NUMBA_CACHE_DIR`、`MPLCONFIGDIR`、`TMPDIR` 指向可写 scratch 子目录。

首期入口脚本只创建 scratch 缓存子目录，然后 `exec "$@"`，默认命令为 `/opt/venv/bin/eacbp --help`；显式命令原样返回退出码。Docker 配置 init 处理子进程，但 init 不替代外部工具终止测试。停止/SIGTERM、STAR/kb/CellBender 子进程退出和半成品隔离需验收。以作业退出码和 summary/report 判断状态，不设置 HTTP HEALTHCHECK。Apptainer exec 不依赖 Docker 入口脚本，Slurm 包装脚本须单独创建同样的缓存目录并传入环境。

## 5. 本地、集群与资源策略

### 本地/CI

Docker 运行 Linux 镜像。下面是 **镜像实现后的预期调用示例**，当前镜像尚不存在；命令为 Linux/WSL shell，宿主目录需事先创建。

```bash
docker run --rm --init \
  --user "$(id -u):$(id -g)" --cpus 4 --memory 16g \
  --mount type=bind,src="$PWD/data",dst=/data,readonly \
  --mount type=bind,src="$PWD/config",dst=/config,readonly \
  --mount type=bind,src="$PWD/refs",dst=/refs,readonly \
  --mount type=bind,src="$PWD/outputs",dst=/outputs \
  --mount type=bind,src="$PWD/scratch",dst=/scratch \
  -e OMP_NUM_THREADS=4 -e MKL_NUM_THREADS=4 \
  -e OPENBLAS_NUM_THREADS=4 -e NUMBA_NUM_THREADS=4 \
  eacbp:analysis-<source-id> \
  /opt/venv/bin/eacbp run \
  --manifest /config/manifest.json --config /config/analysis.json \
  --data /data/study.h5ad --run-dir /outputs/run_001
```

`manifest.json` 使用 `StudyManifest` 字段，`analysis.json` 使用 `RunConfig` 字段，示例见 [ARCHITECTURE_P3.md](ARCHITECTURE_P3.md)。预先创建宿主 outputs/scratch/config/refs/data 目录，但 **不要创建 outputs/run_001**。恢复时使用相同镜像和挂载，将末尾命令换成 `/opt/venv/bin/eacbp resume --run-dir /outputs/run_001`。检查与报告重建分别使用 `inspect --run-dir ...` 和 `report --run-dir ...`。

16 GiB/4 CPU 只是小规模 h5ad 的起测资源，不构成容量保证。FASTQ 可从现有 Slurm 模板的 16 CPU/64 GiB 起测；实际内存按物种索引、细胞/基因数、稀疏矩阵处理和方法测量。特别记录读取、哈希校验、量化与分析峰值内存、耗时和磁盘占用。必要时降低嵌套 BLAS 线程，避免工具线程与库线程相乘。

### Slurm

优先使用 Apptainer，但这是等待集群环境核验的部署决策。构建镜像及转换 SIF 在 CI/专用构建环境完成；计算节点只读取已准备好的 SIF，避免运行时下载。登录节点只做轻量准备与校验，重计算仍由 Slurm 调度。

新模板保留资源申请与日志约定，将宿主 Python 调用改成 `apptainer exec --cleanenv` 后执行 `/opt/venv/bin/eacbp`，旧 FASTQ 脚本则显式使用 `/opt/venv/bin/python`。绑定 data/refs/config/outputs/scratch，并通过 `--env` 或 `APPTAINERENV_*` 传入缓存及线程变量。镜像内线程数与 `SLURM_CPUS_PER_TASK` 对齐。Slurm 日志目录必须在 `sbatch` 前创建。

Apptainer 可运行 OCI 来源镜像，支持显式绑定；其用户和 home 挂载行为与 Docker 不同，因此不得依赖 Docker `USER` 来保证 Apptainer 权限，依赖实际提交用户及目录权限。参见 [Apptainer OCI 兼容说明](https://apptainer.org/docs/user/latest/docker_and_oci.html)。上线需记录 OCI digest 与生成 SIF 的 SHA-256。

`--cleanenv` 不能单独保证隔离宿主 home/cwd 的挂载。按集群版本与管理员配置审查默认绑定，避免宿主 Python 包或配置覆盖镜像路径；规则依据 [Apptainer 绑定目录说明](https://apptainer.org/docs/user/latest/bind_paths_and_mounts.html)。必须在目标计算节点验收，不能只在登录节点验证命令存在。

如果集群不提供 Apptainer，先确认管理员支持的容器运行时；现有宿主 Python/Slurm 路径保留到容器方案验收完成，不自行在集群部署 Docker daemon。

## 6. 断点恢复与可复现边界

当前已支持 `eacbp resume --run-dir ...`，内部复用 Python API 的严格恢复机制。`run_standard_study.py` 的新 UUID 和 FASTQ 脚本的 `--run-id` 不等同于统一 CLI 的恢复接口；旧脚本恢复不作为 P0 承诺。

CLI 已保存规范化配置并复用已导入 artifact，恢复不重新读取外部 h5ad。容器测试应确认原文件移动后仍可恢复、已保存产物或外部依赖资源篡改被拒绝。新建作业与恢复作业明确区分；失败默认交给用户/调度器处理，不能先设置自动重启循环。原始导入或运行元数据未完整落盘的中断单独验证；无法恢复时保留失败目录并另建 run，不承诺任意时刻强杀都可续跑。

恢复验收覆盖同一镜像、同一挂载路径、保存的 manifest/参数、已导入输入与参考哈希。`platform.platform()` 含宿主相关信息，Linux/Windows 历史任务和不同内核节点间不承诺直接复用。跨节点恢复如确有需求，应先定义兼容性规则再调整指纹，不能为了恢复而关闭源码、依赖或输入完整性校验。

指纹还包括当前环境内所有已安装 distributions。因此 `analysis-test` 与 `analysis` 即使科学依赖相同，也不能默认互相恢复；analysis 与 extended/FASTQ 变体亦然。恢复验收在同一个正式 runtime 镜像内完成，测试镜像仅负责独立回归。Docker 与 SIF 即使来自同一 OCI，也需单独验证环境差异。

镜像版本固定也不能保证不同 CPU、线程数和数学库行为逐位相同。数值结果使用明确容差，科学标签、审计门禁和失败语义要求一致。每次运行另存运行命令、镜像/SIF 哈希、参考数据版本、线程/资源配置及宿主信息。

## 7. 实施顺序与交付物

| 阶段 | 工作和拟新增文件 | 退出条件 |
|---|---|---|
| P0：最小分析镜像 | `docker/Dockerfile`、`docker/entrypoint.sh`、`.dockerignore`、`requirements/linux-py311-analysis.lock`、测试锁；新增容器 CI workflow、h5ad/恢复验收脚本 | Linux 锁可安装，安装 wheel 后从仓库外非 root 运行 h5ad；同 runtime 恢复；输出持久化，失败退出非零 |
| P1：扩展、量化与集群 | extended/kb/STAR targets 及各自依赖锁、工具版本清单；`slurm/run_container.sbatch`、运行/存储预检脚本；按需 CellBender target | 高级功能实际计算；真正调用 kb/STAR 完成小型 FASTQ 比对；Slurm 节点运行 SIF；目标存储通过验证 |
| P2：发布与跨节点验证 | 跨节点兼容性结果、部署文档、版本回退步骤、镜像来源/依赖清单及 SBOM | 事务提交前后中断、并发隔离、篡改检测通过；按 digest/SIF 哈希可回退；明确支持和拒绝的迁移组合 |

P0 不依赖跨节点恢复；P1 不把 demo 或 mock 当作真实工具验收。未通过容器验收之前，只能宣称代码支持 CLI 恢复，不能宣称 Docker/Apptainer 恢复已经验证。

实施时同步更新 README 安装入口、SERVER_MIGRATION 的容器契约与环境文件适用说明。首期保留现有 bare-metal 路径，避免改变已使用的作业模板。可先本地构建与验证；镜像发布位置由部署环境确定，不在规划阶段推送到外部仓库。

## 8. 验收矩阵

| 类别 | 必需检查 |
|---|---|
| 镜像基础 | 指定平台构建；`pip check`；standard/fate 导入；包内 CLI 六个子命令 `--help`；从 `/outputs` 执行而非仓库 cwd；兼容 target 另查脚本 |
| 科学与回归 | analysis-test 执行完整现有 pytest；正式 runtime 单独执行 h5ad 小样本，避免测试环境掩盖生产缺依赖 |
| 存储 | `test_storage_regressions.py`、`test_task_transactions.py`、`test_crash_and_concurrency.py` 在目标卷上验证；退出容器后重读 artifact、校验 hash |
| FASTQ | kb 与 STAR 各自真实小数据/匹配索引；多 lane、gzip、缺依赖、无效参考、工具失败均覆盖；明确区分 mocked tests 与实际执行 |
| 高级扩展 | extended 实际执行 PyDESeq2、decoupler、Scrublet、CellTypist、LIANA；模型/网络来源明确；无依赖跳过不得计为覆盖；CellBender 单独验证真实输出和日志 |
| 失败语义 | 输入缺失、权限不足、量化失败、审计失败时非零退出；real 模式不生成 synthetic 替代物 |
| 权限与清理 | 指定 UID/GID；只读输入；可写缓存；终止父进程后无遗留工具进程；清理 scratch 不损坏 artifacts |
| 恢复 | 同 digest 重启复用成功任务；改变数据/参数/依赖拒绝错误复用；提交前后崩溃场景；跨节点按指纹规则处理 |
| 资源 | 记录峰值内存、线程与磁盘；超出调度资源时明确失败，不自动降级成其他方法 |
| 发布 | 来源与依赖清单齐全；按不可变 digest 运行；旧版本与其完整输出一起保留，新版本另开 run 验证 |

容器 pytest 显式指定可写 `--basetemp=/scratch/pytest`，避免项目默认 `.pytest_temp` 写入只读安装路径。科学回归可使用有明确标签的合成 fixture；上线还需固定来源与哈希的实际 h5ad/FASTQ 小样本，记录退出码、任务/审计状态、产物哈希和报告。运行时缺包、静默跳过或自动替换算法不满足验收。

## 9. 实施前待核实事项

1. **构建环境**：提供可工作的 Linux Docker daemon 或 Linux CI runner；本机当前探测不通过，不影响规划交付。
2. **依赖基线**：在 Linux 验证候选 Python 3.11/科学包组合；本次未生成假定可用的锁文件或镜像 digest。
3. **集群能力**：核查 Apptainer 版本、允许的绑定、SIF 存放位置、计算节点架构和输出文件系统语义。
4. **容量与数据**：以计划运行的最大规模评估 RAM、scratch 和持久化配额；确定具有可用权限的小样本、索引和模型。
5. **发布位置**：本地验证可先行；远端镜像仓库、访问权限和保留周期在发布前确定。发布版本保留旧 digest、SIF 和对应锁文件；回退恢复原镜像和对应 run，不强制新旧环境混用。

本次只完成代码和文档检查及规划，没有安装依赖、构建/发布镜像、运行测试、执行真实比对或提交远端 Slurm 作业。
