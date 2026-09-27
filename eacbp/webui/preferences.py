"""User-level directory preferences, independent of code and research data."""
from __future__ import annotations

import json
import os
from pathlib import Path

from pydantic import BaseModel, ConfigDict, Field


class DirectorySettings(BaseModel):
    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True, strict=True)
    workspace: str = Field(min_length=1)
    runs_dir: str = Field(min_length=1)


def absolute_directory(value: str) -> Path:
    path = Path(os.path.expandvars(value)).expanduser()
    if not path.is_absolute():
        raise ValueError("请填写当前运行机器上的完整绝对目录路径")
    return path.resolve()


def default_settings_file() -> Path:
    if os.name == "nt":
        base = Path(os.environ.get("LOCALAPPDATA") or Path.home() / "AppData" / "Local")
    else:
        base = Path(os.environ.get("XDG_CONFIG_HOME") or Path.home() / ".config")
    return base / "eacbp" / "webui.json"


def startup_directories(*, workspace=None, runs_dir=None, settings_file=None):
    """Explicit CLI paths override the corresponding saved preference."""
    path = Path(settings_file or default_settings_file()).expanduser().resolve()
    saved = {}
    if path.exists():
        if path.stat().st_size > 64_000:
            raise ValueError(f"目录设置文件过大：{path}")
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
            if not isinstance(payload, dict) or type(payload.get("schema_version")) is not int or payload["schema_version"] != 1:
                raise ValueError("不支持的目录设置版本")
            saved = DirectorySettings.model_validate({key: value for key, value in payload.items() if key != "schema_version"}).model_dump()
            for key, value in saved.items():
                saved[key] = absolute_directory(value)
        except (ValueError, TypeError) as exc:
            raise ValueError(f"无法读取目录设置 {path}：{exc}") from exc
    data_root = Path(workspace if workspace is not None else saved.get("workspace", Path.cwd())).expanduser().resolve()
    output_root = Path(runs_dir if runs_dir is not None else saved.get("runs_dir", data_root / "outputs" / "runs")).expanduser().resolve()
    if not data_root.is_dir():
        raise ValueError(f"输入数据目录不存在：{data_root}；可用 --workspace 指定现有目录")
    return data_root, output_root, path
