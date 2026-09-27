# P4：统一准入、执行约束与可搬迁研究

本轮保留模块化单体、顺序 DAG 调度、独立科学审计、原子产物提交和严格恢复。重点完善以下六个边界。

| 方向 | 当前实现 | 主要入口 |
|---|---|---|
| 证据准入 | 内置和插件抽取均生成候选证据，由统一入口验证 task/output 身份，沿全部祖先计算 roots 并继承模拟标记 | `eacbp/evidence/provenance.py` |
| 执行约束 | 编排阶段的累计截止时间、取消信号、GPU 禁用；外部命令受进程监督；控制错误不按普通代码错误重试 | `eacbp/execution_context.py`、`orchestrator/execution.py` |
| 能力合同 | 描述符支持 shared/per_target；注册扩展装配后统一展开目标分支；常用方法有严格参数模型 | `capabilities/parameters.py`、`orchestrator/planning.py` |
| 存储读取 | 元数据、流式文件校验和 payload 加载分开；恢复核验不再反序列化矩阵；批量查询减少索引刷新 | `artifact/registry.py`、`artifact/storage.py` |
| 应用服务 | CLI、WebUI worker 和标准研究脚本共用包内公共研究服务，CLI 负责参数与输出 | `application/study_service.py` |
| 研究搬迁 | 注册表 v2/快照 v2 保存相对存储位置；研究包校验文件、计算和审计收据，支持导出/导入后报告重建 | `artifact/portability.py` |

## 使用

原有 `eacbp run/resume/report` 和 `eacbp.cli.run_study` 入口保留。Python 使用新的公共服务：

```python
from eacbp.application.study_service import run_study
from eacbp.execution_context import ExecutionContext

context = ExecutionContext(threads=2)
summary = run_study(manifest=manifest, config=config,
                    data="study.h5ad", run_dir="outputs/run_001",
                    execution_context=context)
# 另一个线程可调用 context.cancel() 请求取消。
```

调用者上下文只能收紧 manifest 的截止时间/GPU 策略。线程设置用于外部进程的常见数值库环境变量；工具自身的线程选项仍由对应参数决定。Python 进程内计算在任务边界检查截止时间和取消，超时返回的产物不会发布；这不是任意原生库调用的硬抢占。外部命令运行中轮询取消与超时，并终止受监督的子进程。尚未增加 WebUI 取消按钮、分布式调度或中间产物保留策略。

恢复复用同样检查取消与截止时间。审计或证据抽取期间发生取消时，不准入本任务的证据，已完成的计算与待审回执保留，后续恢复可以直接重审。

```powershell
eacbp export --run-dir outputs/run_001 --output study.eacbp.zip
eacbp import --bundle study.eacbp.zip --run-dir outputs/imported_001
eacbp report --run-dir outputs/imported_001 --output rebuilt.md
```

研究包包含已注册产物、恢复/事件记录、证据快照和明确的运行说明文件；仅最新 snapshot.json 作为当前可重建证据复核，历史 snapshots 作为保留原字节的历史附件；不复制未提交事务，也不自动收集外部参考库/工具程序。目标文件与导入目录必须不存在。搬迁后可校验并重建报告；继续计算仍要求原配置、方法、源码/依赖指纹和外部资源满足严格恢复条件。元数据哈希与包清单用于检测损坏，不是第三方签名。

应用层导出与导入使用相同的研究包准入检查；缺少完整运行配置、配置版本不支持或 manifest 与最新快照不一致时，导出不会发布 ZIP，导入不会发布目标运行目录。

## 兼容与科学边界

- 旧注册表 v1 和快照 v1 仍可读取。旧绝对路径记录应在原位置导出成新包，再导入新目录；直接搬动旧目录并不保证可读。
- 新快照重新验证真实 roots 和模拟标记。插件不能自报 `source_verified=True`；当前没有独立可信的来源认证回执，明确拒绝此标记。独立审计和所有必需检查仍保留。
- 参数模型校验显式设置，不把新增默认值写入旧合同；未知或拼错参数会提前拒绝。编排上下文有独立声明。方法实现中的科学检查不被参数模型替代。
- `scope="per_target"` 的扩展由框架生成分支与重连产物；`scope="shared"` 的扩展不自动复制。多目标共享扩展引用未展开的目标输出会报错，不能隐式选择某一个目标。
- 配置和源码升级本就会触发严格恢复拒绝；本轮未增加绕过或自动覆盖旧运行的开关。
- 保留文件式索引，未引入 SQLite；已提供批量读取接口，数据库迁移留待实测规模需求。

## 验证

测试覆盖插件来源/模拟传播、同源证据、参数拒绝与多目标接线、进程取消/超时/GPU 策略、旧入口与共享服务、无矩阵加载的恢复校验、相对路径搬迁、包篡改与路径越界、导入后的报告重建。

2026-09-21 在 Windows / Python 3.13 的项目内 `.venv` 安装 `.[dev,standard,fate,advanced-statistics,advanced-qc,communication]` 后，全量回归结果为 **405 passed、1 skipped**（161.11 秒，423 条警告）。唯一跳过项为 `test_cellbender_real_report_keeps_elbo_warnings_beyond_jupyter_css`：当前环境缺少外部 CellBender 冒烟报告。日志与 JUnit 记录位于本地 `.cache/p4/final.log` 和 `.cache/p4/final.xml`，依赖版本保存在 `.cache/p4/environment.txt`。

```powershell
$env:CELLTYPIST_FOLDER = Join-Path (Get-Location) '.cache/celltypist'
$env:NUMBA_CACHE_DIR = Join-Path (Get-Location) '.cache/numba'
.venv/Scripts/python.exe -m pytest tests -q -rs --disable-warnings -p no:cacheprovider --basetemp=.pytest_temp_p4_final
```

新增架构图已通过 Archify 规范检查与四种视口的浏览器检查，并实际审阅深色/浅色截图。Windows 进程树超时与取消已实测；POSIX 进程组分支尚未在 Linux 运行。未运行真实 FASTQ 比对、大规模性能基准或集群部署。
