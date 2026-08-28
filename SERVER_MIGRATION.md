# EACBP 服务器连接、初始化须知与部署指引

本项目严格遵循通用服务器部署模板（`/public/home/qiaoke/remote_python_project_bundle.tar.gz`）及 `Perturb-seq_test_project` 的规范。

---

## 1. 核心运行契约与节点分工 (Execution Contract)

1. **目标平台**：Linux x86-64。
2. **登录节点（Login Node）**：
   - 仅用于准备目标项目专属的 Python 运行环境、同步代码、执行非侵入式健康检查与静态校验；
   - **禁止在登录节点直接运行长时间/高负载的重型计算任务**。
3. **计算节点（Compute Node）**：
   - 由 Slurm 作业调度器统一调度；
   - 计算脚本中必须调用由部署者准备好的 **Python 解释器绝对路径**。
4. **依赖管理原则（External to Bundle）**：
   - 部署包不打包虚拟环境、不含本地操作系统硬编码路径；
   - 服务器端 Python 依赖独立配置并通过国内镜像（清华源/阿里源）加速构建。
5. **安全与数据管理契约**：
   - 严禁将密码、私钥、Token 写入代码、脚本、`.env` 或提交历史中；
   - 原始大文件数据集（raw data）、中间产物（interim/processed）和计算输出（outputs）与源码物理解耦，使用独立目录管理。

---

## 2. 远端服务器与 SSH 配置

在 `~/.ssh/config` 中已配置：
```sshconfig
Host 192.168.201.226
  HostName 192.168.201.226
  User qiaoke
  IdentityFile ~/.ssh/perturb_seq_server
  IdentitiesOnly yes
```

---

## 3. 服务器端部署与非侵入式校验 (Non-mutating Verification)

### 进入工作目录
```bash
ssh 192.168.201.226
cd ~/eacbp_project/current
```

### 环境配置与依赖安装
```bash
PYTHON=/public/home/qiaoke/.local/share/mamba/envs/perturb-seq/bin/python

# 快速安装与本地包可编辑挂载
"$PYTHON" -m pip install "pydantic>=2.5.0" -i https://pypi.tuna.tsinghua.edu.cn/simple
"$PYTHON" -m pip install --no-deps -e .
```

### 登录节点环境校验
```bash
bash scripts/verify_remote.sh
```

---

## 4. Slurm 作业提交与集群运维 (Slurm Cluster Specs)

### 提交计算作业 (CPU 任务模板)
```bash
cd ~/eacbp_project/current
sbatch slurm/run_cpu.sbatch
```

### Slurm 调度参数说明：
- **账号与分区**：`#SBATCH --account=libin`，`#SBATCH --partition=gpu`（集群未配置独立 CPU 分区时，通用任务挂载到 gpu 分区运行）
- **计算资源**：1 节点、1 task、4 cpus-per-task、16GB 内存
- **多线程控制**：
  ```bash
  export CUDA_VISIBLE_DEVICES=""
  export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
  export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
  export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
  ```
- **日志输出路径**：`logs/slurm-%j.out` 与 `logs/slurm-%j.err`

### 常用集群运维与监控命令：
```bash
# 1. 查看作业队列
squeue -u qiaoke

# 2. 查询作业执行资源与状态
sacct -j <JOB_ID> --format=JobID,State,Elapsed,AllocTRES,MaxRSS,ExitCode

# 3. 实时跟踪日志
tail -f logs/slurm-<JOB_ID>.out

# 4. 取消/终止作业
scancel <JOB_ID>
```
