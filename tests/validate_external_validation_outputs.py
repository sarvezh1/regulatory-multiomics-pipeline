#!/usr/bin/env python3
"""Validate retained external-validation evidence."""

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


def rows(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        table = list(csv.reader(handle, delimiter="\t"))
    assert table and len(set(table[0])) == len(table[0]), f"malformed header: {path}"
    assert all(len(row) == len(table[0]) for row in table[1:]), f"malformed row: {path}"
    return [dict(zip(table[0], row)) for row in table[1:]]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default=".")
    args = parser.parse_args()
    root = Path(args.root)
    output = root / "tests/output/external_validation"
    replay1, replay2 = output / "crispr_replay1", output / "crispr_replay2"

    freeze = {int(row["maximum_distance_bp"]): row["sha256"] for row in rows(root / "metadata/external_validation_prediction_freeze.tsv")}
    for distance, expected in freeze.items():
        actual = sha(output / f"full_genome_predictions/predictions/candidate_links.{distance}.tsv")
        assert actual == expected, f"candidate prediction changed after truth comparison: {distance}"

    for first in sorted(replay1.glob("*.tsv")):
        second = replay2 / first.name
        assert second.is_file() and sha(first) == sha(second), f"deterministic replay mismatch: {first.name}"

    denominator = {row["stage"]: int(row["n_records"]) for row in rows(replay1 / "crispr_benchmark_denominators.tsv")}
    assert denominator["total_source_records"] == 10356
    assert denominator["final_scored_records"] == 10288
    assert denominator["experimentally_labelled_positive_records"] == 427
    assert denominator["experimentally_labelled_negative_records"] == 9861

    metrics = rows(replay1 / "crispr_benchmark_metrics.tsv")
    primary = next(row for row in metrics if row["method"] == "candidate_100000" and row["source_study"] == "ALL")
    assert (primary["TP"], primary["FP"], primary["TN"], primary["FN"]) == ("71", "51", "9810", "356")
    assert primary["n_scored"] == "10288"

    pairs = rows(replay1 / "crispr_benchmark_pairs.tsv")
    assert len({row["benchmark_pair_id"] for row in pairs}) == len(pairs)
    assert all(row["coordinate_system"] == "0-based-half-open" for row in pairs)
    failures = Counter(row["classification"] for row in rows(replay1 / "crispr_failure_classes.tsv"))
    assert failures == Counter({"TRUE_POSITIVE": 71, "FALSE_POSITIVE": 51, "TRUE_NEGATIVE": 9810, "FALSE_NEGATIVE": 356})

    ccre = rows(output / "ccre/ccre_concordance_summary.tsv")
    assert {row["ccre_category"] for row in ccre if row["ccre_category"] != "ANY_CCRE"} == {"PLS", "pELS", "dELS", "CA-CTCF", "CA-H3K4me3", "CA-TF", "CA-only", "Low-DNase"}
    statuses = {row["component"]: row["status"] for row in rows(output / "external_validation_status.tsv")}
    assert statuses["crispr_functional_benchmark"] == "COMPLETE"
    assert statuses["optional_differential_validation"] == "NOT_AVAILABLE_INPUT"
    provenance = json.loads((output / "external_validation_provenance.json").read_text(encoding="utf-8"))
    assert provenance["prediction_frozen_before_truth_comparison"] is True
    assert provenance["integration_match_type"] == "REFERENCE_COMPARABLE"
    print("External-validation retained-output validation passed")


if __name__ == "__main__":
    main()
