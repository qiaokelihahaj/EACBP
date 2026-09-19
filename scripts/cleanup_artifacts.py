"""Preview or remove old, unreferenced, inactive EACBP transactions."""
import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from eacbp.artifact.maintenance import cleanup_artifacts, ArtifactMaintenanceError


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact_root", type=Path)
    parser.add_argument("--older-than-days", type=int, default=30)
    parser.add_argument("--apply", action="store_true", help="Delete eligible transactions; default is preview only")
    args = parser.parse_args(argv)
    try:
        report = cleanup_artifacts(args.artifact_root, apply=args.apply, older_than_days=args.older_than_days)
    except (ArtifactMaintenanceError, ValueError) as exc:
        parser.error(str(exc))
    print(json.dumps(report.as_dict(), ensure_ascii=False, indent=2))
    return 1 if report.blocked else 0


if __name__ == "__main__":
    raise SystemExit(main())
