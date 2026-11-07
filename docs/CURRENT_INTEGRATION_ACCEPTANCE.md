# P1 当前集成验收记录

记录日期：2026-09-27。Linux Python 3.12 全量集成测试、通信模块定向测试和最终 wheel 的隔离环境冒烟检查均已通过。Docker Desktop 代理配置修复后，固定摘要镜像已构建，非 root 缓存写入探针通过；合成容器端到端冒烟全部通过（8 个任务，0 个阻断，模拟数据标记已验证）。Slurm 端到端验收未执行。

## 环境和变更范围

Linux 检查在 WSL Ubuntu x86_64、Python 3.12.3 中进行；锁文件由 uv 0.11.8 生成。变更包括：类型化读取继续兼容已有的 `artifact.get/load_payload` 通用读取 API，并保留审计、哈希和完整性检查；H5AD 读写保留元数据；NPZ 读取保留归档中的索引名称。

部署镜像针对 Linux/amd64，包含 `standard`、`fate`、`advanced-statistics`、`advanced-qc`、`communication` 五个科学运行时 extra。第三方运行依赖锁为 97 个固定版本及 PyPI SHA-256，位于 `docker/requirements-linux-py312-scientific.lock`；Python 基础镜像固定为 Python 3.12.14 Bookworm OCI 索引摘要 `sha256:392307d22300de8b5986851a12d9176dfc0fc073e65bf6523ebd7dcbeb23564e`。WSL 环境使用 uv 创建（未包含 pip），安装锁文件及开发测试依赖后，`uv pip check` 通过，共 103 个已安装包。安装日志和 freeze 快照在 `outputs/full-linux-install.log`、`outputs/full-linux-freeze.txt`。

主要验证命令：

```bash
uv pip check --python /tmp/eacbp-linux-full-312/bin/python
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 \
  /tmp/eacbp-linux-full-312/bin/python scripts/verify_ci_test_suite.py integration \
  --junitxml=outputs/full-linux-final.xml
```

## 测试快照

| 检查 | 结果 | 证据 |
|---|---|---|
| Windows Python 3.13 核心测试 | 364 passed，0 skipped，1 deselected；72.45 秒。包含 3 个 runner guard，早于新增的 transaction 测试 | `outputs/core-windows-py313.log`、`outputs/core-windows-py313.xml` |
| Linux Python 3.11 核心测试 | 363 passed，1 个预期的 Windows-only skip，1 deselected；87.23 秒。包含 3 个 runner guard，早于新增的 transaction 测试 | `outputs/core-linux-py311.log`、`outputs/core-linux-py311.xml` |
| Linux Python 3.12 核心测试 | 364 passed，1 个 Windows-only skip，1 个 H5AD 测试 deselected；68.10 秒 | `outputs/core-linux-py312-final.log`、`outputs/core-linux-py312-final.xml` |
| Linux Python 3.12 全量集成测试 | 546 passed，2 个明确 skip，共 548 项；215.66 秒。跳过项为缺少外部 CellBender 真实报告，以及 Windows-only process-assignment 检查 | `outputs/full-linux-final.log`、`outputs/full-linux-final.xml` |
| Linux Python 3.12 `.[dev,standard]` capability lane | 101 passed，0 skipped，1 个 CellRank 用例显式 deselected；88.47 秒。no-skip gate 通过 | `outputs/standard-linux-final.log`、`outputs/standard-linux-install.log` |
| 通信模块定向测试（通信变更后） | 9 passed，0 skipped；19.13 秒 | `outputs/communication-linux-final.log` |
| 最终 wheel 隔离环境冒烟检查 | 通过：CLI、WebUI 资源、验收模板、已安装 Slurm 模板及 deployment CLI；wheel 在最后两项源代码清理后重建 | `outputs/wheel-linux-final.log` |
| 最终 wheel 源码一致性 | wheel 内所有 `eacbp/*.py` 文件与该阶段验收时工作树源码逐项一致（按行尾规范化后比较）；SHA-256：`feb5562156085083f12c6404e029a4a9f45f918d592a374129bd568888c32e12` | `outputs/linux-wheels/eacbp-0.1.0-py3-none-any.whl` |
| runner guard 单独报告 | 3 passed；这 3 项也已包含在 Windows Python 3.13 和 Linux Python 3.11 核心快照中，不应重复计数 | `outputs/ci-guard-windows.xml` |

早期分组前的核心全量基线为 24 failed、474 passed、27 skipped，记录在 `outputs/core-baseline.log`、`outputs/core-baseline.xml`。首次 Linux Python 3.12 全量快照为 540 passed、2 failed、2 skipped，记录在 `outputs/full-linux-baseline.log`、`outputs/full-linux-baseline.xml`；之后修复了 H5AD 跨进程重载断言和 STARsolo 多泳道 slash-key HDF5 写入问题，最终复跑结果如上。Windows 和 Python 3.11 核心快照已包含 3 个 runner guard，但早于新增的 1 个 transaction 测试。

## 合成内存探针

探针输入为 200,000 个细胞、30,000 个基因、200,000 个非零值，伪汇总矩阵为 6 × 30,000。汇总用时 1.57 秒，OS 报告最大 RSS 为 246,584 KiB，采样到的 RSS 增量为 34,004,992 字节。审计输入重建用时 1.67 秒，最大 RSS 为 248,528 KiB，采样增量为 39,636,992 字节。

原始记录在 `outputs/linux-pseudobulk-memory.json`、`outputs/linux-pseudobulk-time.txt`、`outputs/linux-audit-memory.json`、`outputs/linux-audit-time.txt`。这是稀疏合成探针；它没有运行完整审计或模型，也不代表真实密度数据或端到端性能。RSS 采样可能漏过短暂峰值。

## 部署和科学验收边界

初次 `docker build --platform linux/amd64 --file docker/Dockerfile ...` 因 Docker Hub `registry-1.docker.io:443` 超时而无法读取固定基础镜像元数据，记录在 `outputs/docker-build-attempt.txt`。根因是 Docker Desktop 未使用主机已运行的 Clash 代理；Windows 系统代理本身保持关闭。为 Docker Desktop 的当前用户配置 `http://127.0.0.1:7890` 后，固定 Python 3.12.14 Bookworm 摘要可拉取，97 个哈希锁定科学运行依赖安装完成，镜像内 `pip check` 通过。初始镜像运行检查发现 root 创建的 `/scratch` 科学库缓存对 UID `10001` 不可写；Dockerfile 已修正该目录所有权。随后重建的最终标签摘要为 `sha256:f65316b1831456ed008d612c1f239e594c5b026e8ac3ca7d88ba68329af969e7`，默认用户 `10001:10001`，大小 395,346,617 字节。BuildKit 成功重试日志为 `outputs/docker-build-runtime-retry.log`，非 root、`--network none` 缓存写入探针通过。合成端到端容器检查全部通过，最终日志为 `outputs/container-runtime-smoke-final-20260927.log`，详见 [DOCKER_ACCEPTANCE.md](DOCKER_ACCEPTANCE.md)。

没有提供真实生物数据、研究设计或 Slurm 集群，因此没有声称完成真实研究验收、集群存储检查、Apptainer 运行或 Slurm 作业提交。合成测试和方法测试不能替代生物学审查。

## 2026-10-01 公开数据验证后续

PBMC3k 与 Kang 2018 的公开数据、最终结果和隔离容器复跑命令见 [PUBLIC_DATA_VALIDATION.md](PUBLIC_DATA_VALIDATION.md)。最终镜像 `eacbp:py312-cpu` 的 ID 为 `sha256:3cfdb1134215ee4a6b0a524ff61519f95b7f8aa5238dd81ca52ab487af58aa16`，已安装包的 111 个 Python 文件与该阶段验收时工作树逐字节一致。PBMC3k 独立双跑的 PCA、UMAP、邻接图及标签相同，真实运行的恢复、报告和研究包导出/导入均通过；Kang 两个细胞类型的工程、方法和独立审计复跑均通过。Windows 全量 CI 为 587 passed、1 expected skip，Docker/Linux Python 3.12 为 586 passed、2 expected skips，两份 588 项报告的未预期 skip 门禁均通过。数据、缓存和本轮验收记录位于 `F:\EACBP-public-validation`。研究领域审查仍为 `pending_domain_review`；上方旧日期的集成测试快照保留为历史记录。

## 科研价值与效率验收后续

Kang 的供体留一敏感性分析和独立统计参考计算正在单独验收，状态、逐 fold 对照及机器耗时见 [RESEARCH_VALUE_VALIDATION.md](RESEARCH_VALUE_VALIDATION.md)。该参考比较只衡量直接 PyDESeq2 统计计算，不代表完整 GUI/工作流效率基线，也不支持人工时间节省结论；科学解释仍待领域审查。
