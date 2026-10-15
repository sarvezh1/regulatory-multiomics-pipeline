#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument("--manifest", required=True)
p.add_argument("--input-dir", required=True)
p.add_argument("--output", required=True)
a = p.parse_args()

with open(a.manifest, encoding="utf-8", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    rows = list(reader)
    fields = reader.fieldnames
if not rows:
    raise SystemExit("integration manifest is empty")

staged = {}
for item in Path(a.input_dir).iterdir():
    if item.is_file():
        if item.name in staged:
            raise SystemExit(f"duplicate staged basename: {item.name}")
        staged[item.name] = item
for row in rows:
    name = Path(row["path"]).name
    if name not in staged:
        raise SystemExit(f"manifest input was not staged: {name}")
    row["path"] = str(staged[name])

with open(a.output, "w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
