# EACBP 容器化规划

日期：2026-09-16。状态：规划完成，尚未实现或构建镜像。

## 1. 决策摘要

采用 **Linux amd64、CPU 批处理镜像**。本地及 CI 使用 Docker；现有 Slurm 集群优先评估 Apptainer，以同一 OCI 镜像转换为 SIF 执行。一次容器执行对应一个研究作业，结果落到持久化目录。

先交付 h5ad 分析镜像，再增加 FASTQ 工具镜像，最后完善跨容器恢复。现有代码没有 Web 服务、数据库或消息队列，首期无需端口、Compose 服务编排或 Kubernetes。当前科研能力不要求 GPU；Slurm 使用名为 gpu 的分区不等于程序需要 CUDA。

规划默认保留现有计算方法、审计、证据和失败语义；容器化的成功标准是环境与执行可复现，不是改变科学结论。

## 2. 已核实的项目现状

| 代码或文件 | 观察 | 对容器化的影响 |
|---|---|---|
| `pyproject.toml` | Python >=3.10；`standard`、`fate`、`dev` 分组；多数依赖仅设下限 | 必须独立生成 Linux 依赖锁，不能用一次普通 pip 安装作为发布依据 |
| `environment.yml` | Python <3.12，未完整覆盖 standard/fate | 属于旧环境入口，应明确适用范围，不能作为新镜像唯一依赖来源 |
| `requirements-tested-windows-py313.txt` | Windows/Python 3.13 快照，包含 editable Git SSH 安装 | 只作版本参考；禁止直接拿来构建 Linux 发布镜像 |
| `pyproject.toml` 包发现配置 | 排除了 `scripts*`，没有 console scripts | 仅安装 wheel 不会提供现有脚本；首期需要显式复制入口脚本 |
| `scripts/run_standard_study.py` | h5ad 入口，输出使用新 UUID；成功退出 0，失败退出 1 | 可作为首个容器验收入口；目前不能通过 CLI 指定旧 run 恢复 |
| `scripts/run_fastq_to_biology_pipeline.py` | real/demo 显式区分；调用 kb 或 STAR；返回非零失败状态 | FASTQ 镜像必须包含实际工具链，不能用 demo 证明真实比对可用 |
| `eacbp/artifact/` | 使用文件锁、硬链接、原子替换与事务目录 | 持久化卷必须验证这些文件系统语义，不能直接替换成对象存储挂载 |
| `eacbp/orchestrator/checkpoint.py` | 恢复指纹含源代码、依赖和 `platform.platform()` | 镜像相同也不保证跨宿主内核恢复；需要专项验收 |
| `slurm/*.sbatch` | 宿主 Python、项目和日志路径存在绝对路径默认值 | 新增容器模板，保留现有运行方式作为迁移回退 |
| `.github/workflows/tests.yml` | 已配置 Windows/Linux、Python 3.11/3.12/3.13 测试和 wheel 构建 | 可扩展容器测试；配置存在不代表远端测试已成功 |

本地发现 `docker.exe` 与 `wsl.exe`，未检查 Docker daemon/Linux 引擎状态；未连接集群核查 Apptainer、节点权限、存储与配额。工作区已有大量未提交修改，实施时需先选定可追溯源码快照，不能只用当前 HEAD 标识实际镜像内容。

## 3. 镜像边界与依赖策略

### 3.1 发布目标

| 目标 | 内容 | 交付阶段 |
|---|---|---|
| `analysis` | EACBP + standard + fate；支持 h5ad、空间分析和现有本地知识模块 | P0 |
| `analysis-test` | 与 analysis 同一运行依赖，加测试工具和测试文件 | P0，仅用于验证 |
| `fastq-kb` | analysis + 固定版本 kb-python/kallisto/bustools | P1，优先对应当前默认路径 |
| `fastq-star` | analysis + 固定版本 STAR 及 gzip 解压所需命令 | P1，按 STARsolo 用户需求验收 |

不先拆分编排器、审计器和能力模块为多个服务：它们当前共享本地事务注册表与进程内状态，拆分会增加新的协议和一致性工作。

### 3.2 构建约定

1. 候选基线为 Python 3.11 + Debian bookworm slim，目标 `linux/amd64`；这是兼顾现有环境约束的待验证选择。若当前依赖无法解析，应评估调整版本或 Python 3.12，并同步迁移文档，不能静默降级科学方法。
2. 在目标 Linux/Python 环境解析完整传递依赖并保存精确版本和哈希；保留 `harmonypy==0.0.10`、CellRank/AnnData 范围等现有约束。运行与测试使用相同科学计算依赖锁。
3. 基础镜像固定实际可获取的 digest；Python 包构建 wheelhouse，应用安装非 editable wheel。基础系统包与外部工具记录版本、来源、校验和；不使用浮动最新版工具下载。
4. 多阶段构建：builder 负责 wheel/必要编译；runtime 仅保留运行依赖。所需动态库通过安装与导入测试确定，不预设一长串未验证系统包。
5. 首期将三个运行脚本显式复制到 `/opt/eacbp/scripts/`。包与虚拟环境放在 `/opt` 下，不依赖用户 home。后续可增加正式 CLI 包装，避免永久维护脚本复制方案。
6. 禁止 `COPY . .` 无差别打包。应用源码、必要元数据和入口采用白名单 COPY；`.dockerignore` 排除 `.git`、`.agents`、`.venv`、`.artifacts`、所有 pytest/temp 目录、outputs/logs/build、原始数据与本地配置。
7. 镜像记录源码 revision、源码内容哈希、依赖锁哈希、基础镜像 digest 与构建时间。发布镜像从确定的源码快照构建；开发脏工作区构建必须有独立标识。
8. kb 安装后验证实际 `kb`、`kallisto`、`bustools` 可执行且运行时无需下载；STAR 验证版本、动态链接和压缩输入路径。工具版本要与引用索引制作版本一起归档。

Docker 官方建议使用多阶段构建、固定基础镜像 digest 和非 root 用户，参见 [Docker 构建最佳实践](https://docs.docker.com/build/building/best-practices/)。具体科学计算依赖的兼容性仍以本项目 Linux 实测为准。

## 4. 文件系统与执行契约

| 容器路径 | 权限 | 内容与生命周期 |
|---|---|---|
| `/opt/eacbp`、`/opt/venv` | 只读 | 脚本、应用与依赖 |
| `/data` | 只读绑定 | FASTQ/h5ad 原始输入，manifest 内使用容器可见路径 |
| `/refs` | 只读绑定 | kb 索引、t2g、STAR genomeDir、whitelist、可选 GTF |
| `/config` | 只读绑定 | 样本 manifest、markers、terminal states 等运行配置 |
| `/outputs` | 持久化读写绑定 | 完整 run 目录、注册表、事务载荷、journal、报告和比对结果 |
| `/scratch` | 作业级读写绑定 | 可丢弃的缓存和临时文件；使用节点本地磁盘时不提供恢复保证 |

固定容器内路径，允许宿主目录变化。样本 manifest 中的宿主绝对路径不能原样传入；准备容器版本 manifest 并保留原文件。所有输入路径在运行前检查可读，输出路径检查可写。

**必须保存整个 run 目录。** `_transactions` 中包含已经提交的真实 artifact，不可视为普通缓存批量删除。索引、载荷与 journal 应在同一个支持所需原子操作的文件系统中；快照或备份在作业停止/一致性点进行。

存储上线前在实际挂载点运行锁、硬链接、原子替换和多进程冲突测试。NFS/并行文件系统、Windows 绑定目录均不能只凭名称假设兼容；不满足时先使用经验证的 Linux 存储。每个研究使用独立 run 目录，不把当前注册表当作分布式数据库。

Docker 使用非 root 用户并支持宿主 UID/GID 对齐；输出目录由部署者预先准备。默认工作目录 `/outputs`，所有脚本仍显式接收 `--output-dir`。建议 `PYTHONUNBUFFERED=1`、`PYTHONDONTWRITEBYTECODE=1`、`PYTHONNOUSERSITE=1`，并把 `HOME`、`XDG_CACHE_HOME`、`NUMBA_CACHE_DIR`、`MPLCONFIGDIR`、`TMPDIR` 指向可写 scratch 子目录。

首期采用直接执行命令、原样返回退出码的入口；Docker 配置 init 处理子进程。停止/SIGTERM、STAR/kb 子进程退出和半成品隔离需验收。以作业退出码和 summary/report 判断状态，不设置无意义的 HTTP HEALTHCHECK。

## 5. 本地、集群与资源策略

### 本地/CI

Docker 运行 Linux 镜像。下面是 **镜像实现后的预期调用示例**，当前镜像尚不存在；命令为 Linux/WSL shell，宿主目录需事先创建。

```bash
docker run --rm --init \
  --user "$(id -u):$(id -g)" --cpus 4 --memory 16g \
  --mount type=bind,src="$PWD/data",dst=/data,readonly \
  --mount type=bind,src="$PWD/outputs",dst=/outputs \
  --mount type=bind,src="$PWD/scratch",dst=/scratch \
  -e OMP_NUM_THREADS=4 -e MKL_NUM_THREADS=4 \
  -e OPENBLAS_NUM_THREADS=4 -e NUMBA_NUM_THREADS=4 \
  eacbp:analysis-<source-id> \
  /opt/venv/bin/python /opt/eacbp/scripts/run_standard_study.py \
  --data /data/study.h5ad --species homo_sapiens --tissue kidney \
  --study-id kidney_study --output-dir /outputs/runs
```

16 GiB/4 CPU 只是小规模 h5ad 的起测资源，不构成容量保证。FASTQ 可从现有 Slurm 模板的 16 CPU/64 GiB 起测；实际内存按物种索引、细胞/基因数、稀疏矩阵处理和方法测量。特别记录读取、哈希校验、量化与分析峰值内存、耗时和磁盘占用。必要时降低嵌套 BLAS 线程，避免工具线程与库线程相乘。

### Slurm

优先使用 Apptainer，但这是等待集群环境核验的部署决策。构建镜像及转换 SIF 在 CI/专用构建环境完成；计算节点只读取已准备好的 SIF，避免运行时下载。登录节点只做轻量准备与校验，重计算仍由 Slurm 调度。

新模板保留资源申请与日志约定，将宿主 Python 调用改成 `apptainer exec --cleanenv` 后执行 `/opt/venv/bin/python`，显式绑定 data/refs/config/outputs/scratch，并通过 `--env` 或 `APPTAINERENV_*` 传入线程变量。镜像内线程数与 `SLURM_CPUS_PER_TASK` 对齐。Slurm 日志目录必须在 `sbatch` 前创建。

Apptainer 可运行 OCI 来源镜像，支持显式绑定；其用户和 home 挂载行为与 Docker 不同，因此不得依赖 Docker `USER` 来保证 Apptainer 权限，依赖实际提交用户及目录权限。参见 [Apptainer OCI 兼容说明](https://apptainer.org/docs/user/latest/docker_and_oci.html)。上线需记录 OCI digest 与生成 SIF 的 SHA-256。

如果集群不提供 Apptainer，先确认管理员支持的容器运行时；现有宿主 Python/Slurm 路径保留到容器方案验收完成，不自行在集群部署 Docker daemon。

## 6. 断点恢复与可复现边界

当前恢复能力在 Python API 的 `current_state={..., "resume": True}`，并非已有通用 CLI 功能。`run_standard_study.py` 总是新建 UUID；FASTQ 的 `--run-id` 也不能等同于恢复开关。

后续需给入口增加稳定 run 标识、显式 `--resume` 和参数匹配检查，复用已有输入 artifact，避免重复注册不可变版本。新建作业与恢复作业必须明确区分；失败默认交给用户/调度器处理，不能先设置自动重启循环。

恢复验收覆盖同一镜像、同一挂载路径、同一 manifest/参数、相同输入与参考哈希。`platform.platform()` 含宿主相关信息，Linux/Windows 历史任务和不同内核节点间不承诺直接复用。跨节点恢复如确有需求，应先定义兼容性规则再调整指纹，不能为了恢复而关闭源码、依赖或输入完整性校验。

镜像版本固定也不能保证不同 CPU、线程数和数学库行为逐位相同。数值结果使用明确容差，科学标签、审计门禁和失败语义要求一致。每次运行另存运行命令、镜像/SIF 哈希、参考数据版本、线程/资源配置及宿主信息。

## 7. 实施顺序与交付物

| 阶段 | 工作和拟新增文件 | 退出条件 |
|---|---|---|
| P0：最小分析镜像 | `docker/Dockerfile`、`.dockerignore`、`requirements/linux-py311-analysis.lock`、测试锁；复制现有脚本；新增容器 CI workflow | Linux 锁可安装，非 root 镜像运行真实 h5ad 小样本，输出离开容器仍可读，失败退出非零 |
| P1：量化与集群 | Dockerfile 增加 kb/STAR targets、独立工具版本清单；`slurm/run_container.sbatch`、运行/存储预检脚本 | 真正调用 kb/STAR 完成小型 FASTQ 比对；Slurm 计算节点运行 SIF；实际挂载点存储测试通过 |
| P2：恢复和发布 | CLI run/resume、容器恢复测试、部署文档与版本回退步骤、镜像来源/依赖清单 | 杀进程后恢复、并发拒绝/隔离、篡改检测通过；按 digest/SIF 哈希可回退 |

P0 不依赖跨节点恢复；P1 不把 demo 或 mock 当作真实工具验收。P2 完成之前文档必须明确“容器可以重跑，但现有 CLI 尚不支持复用旧作业恢复”。

实施时同步更新 README 安装入口、SERVER_MIGRATION 的容器契约与环境文件适用说明。首期保留现有 bare-metal 路径，避免改变已使用的作业模板。可先本地构建与验证；镜像发布位置由部署环境确定，不在规划阶段推送到外部仓库。

## 8. 验收矩阵

| 类别 | 必需检查 |
|---|---|
| 镜像基础 | 指定平台构建；`pip check`；standard/fate 导入；三个入口 `--help`；从 `/outputs` 执行而非仓库 cwd |
| 科学与回归 | analysis-test 执行完整现有 pytest；正式 runtime 单独执行 h5ad 小样本，避免测试环境掩盖生产缺依赖 |
| 存储 | `test_storage_regressions.py`、`test_task_transactions.py`、`test_crash_and_concurrency.py` 在目标卷上验证；退出容器后重读 artifact、校验 hash |
| FASTQ | kb 与 STAR 各自真实小数据/匹配索引；多 lane、gzip、缺依赖、无效参考、工具失败均覆盖；明确区分 mocked tests 与实际执行 |
| 失败语义 | 输入缺失、权限不足、量化失败、审计失败时非零退出；real 模式不生成 synthetic 替代物 |
| 权限与清理 | 指定 UID/GID；只读输入；可写缓存；终止父进程后无遗留工具进程；清理 scratch 不损坏 artifacts |
| 恢复 | 同 digest 重启复用成功任务；改变数据/参数/依赖拒绝错误复用；提交前后崩溃场景；跨节点按指纹规则处理 |
| 资源 | 记录峰值内存、线程与磁盘；超出调度资源时明确失败，不自动降级成其他方法 |
| 发布 | 来源与依赖清单齐全；按不可变 digest 运行；旧版本与其完整输出一起保留，新版本另开 run 验证 |

本次只完成代码和文档检查及规划，没有安装依赖、构建/发布镜像、运行测试、执行真实比对或提交远端 Slurm 作业。
