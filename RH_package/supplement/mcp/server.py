"""Standalone replacements for the certificate's two Rethlas logging calls.

These functions only append local progress records. They do not supply any
mathematical inputs, perform verification, or start a network service.
"""
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def memory_append(problem_id, channel, record):
    directory = ROOT / "run_logs" / problem_id
    directory.mkdir(parents=True, exist_ok=True)
    result = {"timestamp_utc": datetime.now(timezone.utc).isoformat(),
              "channel": channel, "record": record}
    with (directory / (channel + ".jsonl")).open("a", encoding="utf-8") as stream:
        stream.write(json.dumps(result) + "\n")
    return result


def branch_update(problem_id, branch_id, state):
    return memory_append(problem_id, "branch_states",
                         {"branch_id": branch_id, "state": state})
