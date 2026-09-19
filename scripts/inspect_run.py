"""Summarize one EACBP invocation's execution event log."""
import argparse
from collections import Counter
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from eacbp.orchestrator.events import read_run_events


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("events", type=Path, help="JSONL path returned in run_study()['event_log']")
    args = parser.parse_args(argv)
    events = read_run_events(args.events)
    counts = Counter(event.kind for event in events)
    summary = {
        "run_id": events[0].run_id if events else None,
        "study_id": events[0].study_id if events else None,
        "last_event": events[-1].kind if events else None,
        "elapsed_seconds": events[-1].elapsed_seconds if events else None,
        "event_counts": dict(counts),
        "issues": [event.model_dump(mode="json") for event in events
                   if event.kind in {"task_failed", "task_blocked", "audit_rejected", "run_failed"}],
    }
    print(json.dumps(summary, indent=2, ensure_ascii=False, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
