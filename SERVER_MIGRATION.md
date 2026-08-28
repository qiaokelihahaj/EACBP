# EACBP 服务器连接与部署指引

本项目参考了 `Perturb-seq_test_project` 的服务器配置与迁移流程规范。

## 1. 远端服务器与 SSH 配置

在 `~/.ssh/config` 中已配置 SSH 别名与专有私钥：

```sshconfig
Host 192.168.201.226
  HostName 192.168.201.226
  User qiaoke
  IdentityFile ~/.ssh/perturb_seq_server
  IdentitiesOnly yes
```

### 免密连接测试
```powershell
ssh 192.168.201.226 'uname -a; /usr/local/bin/micromamba --version'
```

---

## 2. 代码打包与同步

### 方式 A：通过 Git 仓库同步（推荐）
在服务器目标目录下直接克隆或拉取：
```bash
ssh 192.168.201.226
mkdir -p ~/eacbp_project && cd ~/eacbp_project
git clone https://github.com/qiaokelihahaj/EACBP.git current
cd current
```

### 方式 B：本地打包并 SCP 同步
在 Windows 本地执行：
```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\prepare_server_bundle.ps1
scp -r .\dist\server_bundle 192.168.201.226:~/eacbp_project/current
```

---

## 3. 服务器端环境创建与验证

### 快速接入（使用已验证的共享环境）：
服务器上已配置并安装好依赖的 Python 环境：
```bash
PYTHON=/public/home/qiaoke/.local/share/mamba/envs/perturb-seq/bin/python

cd ~/eacbp_project/current
"$PYTHON" -m pip install "pydantic>=2.5.0" -i https://pypi.tuna.tsinghua.edu.cn/simple
"$PYTHON" -m pip install --no-deps -e .

# 验证环境与运行完整测试套件（92 passed）
"$PYTHON" scripts/verify_environment.py
"$PYTHON" -m pytest -q
```
