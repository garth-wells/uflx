"""Check that test coverage is high."""

import json

with open(".coverage.json") as f:
    data = json.load(f)

too_low = [
    (filename, filedata["summary"]["percent_covered_display"])
    for filename, filedata in data["files"].items()
    if filedata["summary"]["percent_covered"] < 70
]
if len(too_low) > 0:
    raise RuntimeError(
        "Coverage is too low for these files:\n"
        + "\n".join(f"  - {filename} ({percent}%)" for filename, percent in too_low)
    )

if data["totals"]["percent_covered"] < 85:
    raise RuntimeError("Overall coverage is too low")
