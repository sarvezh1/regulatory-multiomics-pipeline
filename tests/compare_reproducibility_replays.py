#!/usr/bin/env python3
import argparse
import hashlib
import json
from pathlib import Path

SCIENTIFIC_PATTERNS = (
    "*.analysis_status.tsv",
    "*.raw_counts.tsv",
    "wgbs.cpg.long.tsv",
    "promoter_evidence.tsv",
    "gene_summary.tsv",
    "*_promoter.tsv",
    "*distal*.tsv",
    "*sensitivity*.tsv",
    "*qc_metrics.tsv",
)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def inventory(root):
    found = {}
    for pattern in SCIENTIFIC_PATTERNS:
        for path in root.rglob(pattern):
            if path.is_file():
                found[path.relative_to(root).as_posix()] = digest(path)
    return found


def normalized_run(path):
    data = json.loads(path.read_text(encoding="utf-8"))
    for key in ("start_timestamp_utc", "end_timestamp_utc", "report_generation_seconds"):
        data.pop(key, None)
    data.get("manifest_checksums", {}).pop("output_manifest.tsv", None)
    return data


parser = argparse.ArgumentParser()
parser.add_argument("first_results")
parser.add_argument("second_results")
parser.add_argument("first_report")
parser.add_argument("second_report")
args = parser.parse_args()

first = inventory(Path(args.first_results))
second = inventory(Path(args.second_results))
if not first or first != second:
    missing = sorted(set(first) ^ set(second))
    changed = sorted(key for key in set(first) & set(second) if first[key] != second[key])
    raise SystemExit(f"MISMATCH scientific outputs: missing={missing}, changed={changed}")

report_a = Path(args.first_report)
report_b = Path(args.second_report)
stable = ["source_manifest.tsv", "input_manifest.tsv", "reference_manifest.tsv", "tool_manifest.tsv", "container_manifest.tsv", "parameter_manifest.json", "status_manifest.tsv", "warning_manifest.tsv", "contract_versions.tsv"]
changed = [name for name in stable if digest(report_a / name) != digest(report_b / name)]
if changed:
    raise SystemExit(f"MISMATCH normalized provenance: {changed}")
if normalized_run(report_a / "run_manifest.json") != normalized_run(report_b / "run_manifest.json"):
    raise SystemExit("MISMATCH normalized run manifests")
print(f"IDENTICAL scientific content: {len(first)} files")
print("SCIENTIFICALLY_EQUIVALENT_EXPECTED_RUNTIME_VARIATION: operational reports and runtime fields excluded")
