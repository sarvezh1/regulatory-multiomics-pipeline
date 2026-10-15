#!/usr/bin/env python3
import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

PRIVATE = [re.compile(r"[A-Za-z]:[\\/]"), re.compile(r"(?:^|[\s\"'])/home/[^/\s]+/"), re.compile(r"(?i)(password|secret|token|api[_-]?key)\s*[:=]")]

def sha(path):
    h=hashlib.sha256(); h.update(path.read_bytes()); return h.hexdigest()

def rows(path):
    with path.open(encoding="utf-8", newline="") as handle: table=list(csv.reader(handle, delimiter="\t"))
    if not table or not table[0] or len(set(table[0])) != len(table[0]): raise ValueError(f"reporting: malformed TSV header: {path.name}")
    for line_number,row in enumerate(table[1:],2):
        if len(row)!=len(table[0]): raise ValueError(f"reporting: malformed TSV row {line_number}: {path.name}")
    return [dict(zip(table[0],row)) for row in table[1:]]

p=argparse.ArgumentParser(); p.add_argument("report_dir"); a=p.parse_args(); root=Path(a.report_dir)
required=["run_manifest.json","input_manifest.tsv","reference_manifest.tsv","tool_manifest.tsv","container_manifest.tsv","parameter_manifest.json","source_manifest.tsv","output_manifest.tsv","status_manifest.tsv","warning_manifest.tsv","report_summary.tsv","run_report.html","contract_versions.tsv"]
for name in required:
    path=root/name
    if not path.is_file() or path.stat().st_size==0: raise ValueError(f"reporting: required artifact missing or empty: {name}")
run=json.loads((root/"run_manifest.json").read_text(encoding="utf-8"))
for name,expected in run["manifest_checksums"].items():
    if sha(root/name)!=expected: raise ValueError(f"reporting: checksum mismatch for {name}")
statuses=rows(root/"status_manifest.tsv")
if not statuses: raise ValueError("reporting: empty status manifest")
summary={r["metric"]:r for r in rows(root/"report_summary.tsv")}
if int(summary["component_count"]["value"])!=len(statuses): raise ValueError("reporting: report component count is not traceable to status_manifest.tsv")
outputs=rows(root/"output_manifest.tsv")
output_ids=[(row["component"],row["relative_path"]) for row in outputs]
if len(output_ids)!=len(set(output_ids)): raise ValueError("reporting: duplicate output identifier")
if int(summary["output_count"]["value"])!=len(outputs): raise ValueError("reporting: report output count is not traceable to output_manifest.tsv")
inputs=rows(root/"input_manifest.tsv")
warnings=rows(root/"warning_manifest.tsv")
if int(summary["input_record_count"]["value"])!=len(inputs): raise ValueError("reporting: report input count is not traceable to input_manifest.tsv")
if int(summary["warning_count"]["value"])!=len(warnings): raise ValueError("reporting: report warning count is not traceable to warning_manifest.tsv")
for modality in ("rna","atac","wgbs","h3k27ac","integration"):
    expected_status=sum(1 for row in statuses if row["component"].lower().startswith(modality))
    expected_outputs=sum(1 for row in outputs if row["component"].lower()==modality)
    if int(summary[f"{modality}_status_component_count"]["value"])!=expected_status: raise ValueError(f"reporting: {modality} status summary is not traceable")
    if int(summary[f"{modality}_output_artifact_count"]["value"])!=expected_outputs: raise ValueError(f"reporting: {modality} output summary is not traceable")
for path in root.rglob("*"):
    if not path.is_file(): continue
    text=path.read_text(encoding="utf-8",errors="ignore")
    if any(pattern.search(text) for pattern in PRIVATE): raise ValueError(f"reporting: private path or credential-like text detected in {path.name}")
print(f"Run report validation passed: {len(statuses)} statuses, {len(outputs)} outputs")
