# Docker 镜像验收

记录日期：2026-09-27。

## 代理根因与恢复

Docker Desktop 的 Linux 引擎直接访问 Docker Hub 时，连接
`registry-1.docker.io:443` 超时。Windows 系统代理处于关闭状态；主机上已有
Clash HTTP 代理监听 `127.0.0.1:7890`，但 Docker Desktop 没有配置使用它。

本次只调整当前用户的 Docker Desktop 代理：停止 Docker Desktop 前确认没有运行中的容器，备份
`%APPDATA%\Docker\settings-store.json` 到
`%APPDATA%\Docker\settings-store.before-eacbp-proxy-20260927.json`，然后将
`ProxyHTTPMode` 设为 `manual`，并将 `OverrideProxyHTTP`、
`OverrideProxyHTTPS` 设为 `http://127.0.0.1:7890`。随后重启 Docker Desktop。
没有修改 Windows 系统代理、WinHTTP 代理或应用代码。拉取镜像期间需要保持
Clash 代理在该地址监听。

如需恢复原 Docker Desktop 设置，先退出 Docker Desktop，另存当前
`settings-store.json`，再从上述备份恢复该文件并重新启动 Docker Desktop。

## 镜像构建

镜像基于固定的 Python 3.12.14 Bookworm OCI 摘要
`sha256:392307d22300de8b5986851a12d9176dfc0fc073e65bf6523ebd7dcbeb23564e`，
安装 `docker/requirements-linux-py312-scientific.lock` 中 97 个带哈希的精确运行依赖，
并用单独的构建锁生成 wheel。镜像标签为 `eacbp:py312-cpu`。

代理修复后的首次 BuildKit 构建成功读取固定基础镜像、安装依赖、通过镜像内
`pip check`，并导出 `eacbp:py312-cpu`。容器运行检查随后发现构建时以 root 导入科学库会在
`/scratch` 下创建 root 所有的缓存，UID `10001` 的运行用户无法写入。Dockerfile 已在构建
导入后增加 `/scratch` 所有权修正。第一次修复后重建在导出阶段遇到 BuildKit snapshot 错误，
记录于 `outputs/docker-build-runtime-fixed.log`；随后重试成功，日志
`outputs/docker-build-runtime-retry.log` 到达 `#20 DONE` 并重新生成 `eacbp:py312-cpu`。此标签
以修复后的构建为准，初始镜像不作为最终验收结果。最终镜像摘要为
`sha256:f65316b1831456ed008d612c1f239e594c5b026e8ac3ca7d88ba68329af969e7`，大小为
395,346,617 字节，默认用户为 `10001:10001`。

## 运行验收

最终镜像通过两个非 root 运行探针：默认 scratch 配置可写；只读根文件系统下的
`/scratch` 临时文件系统和 `/outputs` 输出挂载可写，临时目录为 `/scratch/tmp`。
全部容器禁用网络，输入以只读方式挂载。合成端到端检查通过，退出码为 0：

```bash
python scripts/verify_container_runtime.py \
  --image eacbp:py312-cpu \
  --work-dir outputs/container-runtime-smoke-new
```

复跑时必须指定新的或空的 `--work-dir`；已有证据目录会被拒绝覆盖。

检查确认合成输入文件的 SHA-256 在执行前后不变；运行得到 8 个成功任务、0 个阻塞任务，且
结果标记 `is_simulated=true`。`run`、`resume`、`report`、`export`、`import` 和恢复后的报告重建
全部通过。合成 run ID 为 `f96004a91c234e3f8bb18b7c46b879b9`，输入 SHA-256 为
`4d9bff2052aa81f9313e4a934a020851bb06c0781e30a00ebfa76fa1ca13f02d`。日志和完整证据目录分别为
`outputs/container-runtime-smoke-final-20260927.log` 和
`outputs/container-runtime-smoke-final-20260927/`。

这是容器打包、权限和合成 CLI 生命周期的工程验收，不代表真实生物数据研究验收。Slurm、
Apptainer、集群挂载策略和计算节点存储仍需在目标集群验证。

## 2026-10-01 最终镜像与公开数据验收

本轮最终标签仍为 `eacbp:py312-cpu`，镜像 ID 更新为
`sha256:3cfdb1134215ee4a6b0a524ff61519f95b7f8aa5238dd81ca52ab487af58aa16`，
大小为 395,395,394 字节，Python 为 3.12.14，默认用户为 `10001:10001`。
已安装包的 111 个 Python 源文件与该阶段验收时工作树逐字节一致；revision 标签含
`public-seeded-validation-dirty`，明确记录未提交修改。上述 2026-09-27 镜像和日志保留为历史证据。

Dockerfile 将锁定科学依赖安装移到应用 wheel COPY 之前，revision ARG/LABEL 移到末尾，
并将 pip 下载超时设为 120 秒。最终构建成功，随后只改变 revision 的复建命中依赖和应用
安装缓存。构建日志位于 `F:\EACBP-public-validation\runs\docker-cache-final-build.log`，
缓存探针日志为同目录下的 `docker-revision-cache-probe.log`。

最终镜像再次通过非 root、默认 scratch、只读 root filesystem 及合成 CLI 全生命周期检查，
日志为 `F:\EACBP-public-validation\runs\container-seeded-final.log`，退出码为 0。
合成 run ID 为 `4aaef50f4bfe49d69865afcc7132d47d`，仍标记 `is_simulated=true`。

公开 PBMC3k 与 Kang 2018 复跑使用镜像内已安装包、只读输入、禁用网络和 F: 输出/cache
挂载。PBMC 两个独立容器的 PCA、UMAP、图和标签一致，真实研究包恢复及导出/导入通过；
Kang 的两个供体配对分析和独立审计通过。最终 Docker/Linux 全量 CI 为 586 passed、
2 expected skips，未预期 skip 门禁通过。固定输入、证据路径和完整复跑命令见
[PUBLIC_DATA_VALIDATION.md](PUBLIC_DATA_VALIDATION.md)。研究解释仍需领域审查。
