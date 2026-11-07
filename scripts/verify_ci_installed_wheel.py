"""Smoke-test a built wheel's installed CLI and bundled package resources."""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
import tempfile
import venv
from pathlib import Path


def _run(command: list[str], *, cwd: Path) -> None:
    subprocess.run(command, cwd=cwd, check=True)


def _write_pin_only_constraints(lock_path: Path, destination: Path) -> None:
    pins: list[str] = []
    for original_line in lock_path.read_text(encoding="utf-8").splitlines():
        if not original_line.strip() or original_line.lstrip().startswith("#"):
            continue
        if original_line[:1].isspace():
            continue
        requirement = original_line.strip()
        if requirement.endswith("\\"):
            requirement = requirement[:-1].rstrip()
        if not re.match(r"^[A-Za-z0-9_.-]+==[^\s;]+(?:\s*;.*)?$", requirement):
            raise ValueError(f"unsupported line in hashed lock: {original_line}")
        pins.append(requirement)
    if not pins:
        raise ValueError(f"hashed lock contains no exact package pins: {lock_path}")
    destination.write_text("\n".join(pins) + "\n", encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("wheel", type=Path)
    parser.add_argument(
        "--lock",
        type=Path,
        help="hashed lock whose package pins constrain the independent wheel install",
    )
    args = parser.parse_args(argv)
    wheel = args.wheel.resolve()
    if not wheel.is_file() or wheel.suffix != ".whl":
        parser.error(f"wheel does not exist: {wheel}")
    lock = args.lock.resolve() if args.lock else None
    if lock is not None and not lock.is_file():
        parser.error(f"lock file does not exist: {lock}")

    with tempfile.TemporaryDirectory(prefix="eacbp-installed-wheel-") as temp_dir:
        smoke_root = Path(temp_dir)
        environment = smoke_root / "venv"
        venv.EnvBuilder(with_pip=True).create(environment)
        constraint = smoke_root / "pin-only-constraints.txt"
        if lock is not None:
            _write_pin_only_constraints(lock, constraint)
        if sys.platform == "win32":
            scripts = environment / "Scripts"
            python = scripts / "python.exe"
            cli = scripts / "eacbp.exe"
        else:
            scripts = environment / "bin"
            python = scripts / "python"
            cli = scripts / "eacbp"

        install = [str(python), "-m", "pip", "install"]
        if lock is not None:
            install.extend(["--constraint", str(constraint)])
        install.append(str(wheel))
        _run(install, cwd=smoke_root)
        _run([str(python), "-m", "pip", "check"], cwd=smoke_root)
        _run([str(cli), "--help"], cwd=smoke_root)
        resource_check = (
            "from importlib.resources import files; import sys; from pathlib import Path; "
            "import eacbp; "
            "assert Path(eacbp.__file__).resolve().is_relative_to(Path(sys.prefix).resolve()), eacbp.__file__; "
            "assert files('eacbp.webui').joinpath('static', 'index.html').is_file(); "
            "assert files('eacbp.acceptance').joinpath('fixtures', 'research_manifest.template.json').is_file(); "
            "from eacbp.deployment.slurm import default_batch_script; "
            "batch = default_batch_script(); "
            "assert batch.is_file() and batch.resolve().is_relative_to(Path(sys.prefix).resolve()), batch"
        )
        _run([str(python), "-c", resource_check], cwd=smoke_root)
        _run([str(python), "-m", "eacbp.deployment", "--help"], cwd=smoke_root)

    print("Installed wheel CLI and bundled resources passed the smoke check.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
