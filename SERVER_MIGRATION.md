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
# 进入工作目录
mkdir -p ~/eacbp_project && cd ~/eacbp_project
git clone git@github.com:qiaokelihahaj/EACBP.git current
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

登录服务器后：
```bash
ssh 192.168.201.226
cd ~/eacbp_project/current

# 使用 micromamba 创建虚拟环境
/usr/local/bin/micromamba create -n eacbp python=3.11 -c conda-forge -y
/usr/local/bin/micromamba run -n eacbp pip install -e ".[bio,dev]"

# 验证环境与运行测试套件
/usr/local/bin/micromamba run -n eacbp python scripts/verify_environment.py
/usr/local/bin/micromamba run -n eacbp pytest -q
```
