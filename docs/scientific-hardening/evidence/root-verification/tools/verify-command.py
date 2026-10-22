"""Controller-only raw verification capture; not shipped as a product tool."""
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time
from datetime import datetime, timezone

label, *command = sys.argv[1:]
if not re.fullmatch(r"[a-z0-9-]+", label) or not command:
    raise SystemExit("usage: verify-command.py LABEL COMMAND ...")
if command[0] != "rtk":
    raise SystemExit("verification commands must start with rtk")
out = Path(__file__).parent / "root-verification"
out.mkdir(exist_ok=True)
commit = subprocess.check_output(["rtk", "proxy", "git", "rev-parse", "HEAD"], text=True).strip()
status = subprocess.check_output(["rtk", "proxy", "git", "status", "--porcelain"], text=True).strip()
start = time.monotonic()
result = subprocess.run(command, capture_output=True, env=os.environ.copy())
for suffix, data in (("stdout", result.stdout), ("stderr", result.stderr)):
    (out / f"{label}.{suffix}").write_bytes(data)
text = (result.stdout + result.stderr).decode("utf-8", errors="replace")
counts = [list(map(int, m)) for m in re.findall(r"test result: .*? (\d+) passed; (\d+) failed; (\d+) ignored;", text)]
record = dict(label=label, command=command, returncode=result.returncode,
              source_commit=commit, source_status=status,
              elapsed_seconds=time.monotonic()-start,
              created_utc=datetime.now(timezone.utc).isoformat(),
              rust_test_counts=counts,
              stdout_file=f"{label}.stdout", stderr_file=f"{label}.stderr")
(out / f"{label}.json").write_text(json.dumps(record, indent=2) + "\n")
print(json.dumps(record, indent=2))
print("\n".join(text.splitlines()[-18:]))
raise SystemExit(result.returncode)
