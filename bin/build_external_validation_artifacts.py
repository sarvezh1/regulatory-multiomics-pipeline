#!/usr/bin/env python3
"""Build deterministic external-validation status, provenance, plot sources, and SVG summaries."""

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path


def read_tsv(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, fields, rows):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader(); writer.writerows(rows)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def bar_svg(path, title, rows, label_key, value_key, maximum=None):
    maximum = maximum or max(float(row[value_key]) for row in rows) or 1
    width, left, top, bar_h, gap = 900, 260, 55, 24, 16
    height = top + len(rows) * (bar_h + gap) + 45
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
             '<rect width="100%" height="100%" fill="#ffffff"/>',
             f'<text x="20" y="30" font-family="sans-serif" font-size="20" fill="#17232b">{title}</text>']
    for i, row in enumerate(rows):
        y = top + i * (bar_h + gap)
        value = float(row[value_key]); length = max(0, (width - left - 90) * value / maximum)
        label = str(row[label_key]).replace("&", "&amp;")
        parts.extend([f'<text x="20" y="{y + 18}" font-family="sans-serif" font-size="13" fill="#17232b">{label}</text>',
                      f'<rect x="{left}" y="{y}" width="{length:.2f}" height="{bar_h}" fill="#5b4bb7"/>',
                      f'<text x="{left + length + 8:.2f}" y="{y + 18}" font-family="sans-serif" font-size="12" fill="#17232b">{value:.3f}</text>'])
    parts.append('</svg>')
    Path(path).write_text("\n".join(parts) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True)
    parser.add_argument("--sources", required=True)
    parser.add_argument("--candidate-freeze", required=True)
    parser.add_argument("--ccre", required=True)
    parser.add_argument("--crispr", required=True)
    parser.add_argument("--plausibility", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()
    root, out = Path(args.root).resolve(), Path(args.outdir)

    statuses = [
        {"component": "external_source_audit", "status": "COMPLETE", "reason": "Pinned sources, checksums, roles, and independence classes recorded"},
        {"component": "ccre_concordance", "status": "COMPLETE", "reason": "Current ENCODE K562 annotation concordance completed"},
        {"component": "crispr_functional_benchmark", "status": "COMPLETE", "reason": "Frozen predictions scored against audited positive and negative CRISPR pairs"},
        {"component": "biological_plausibility", "status": "COMPLETE", "reason": "Reference-comparable expression-stratum analysis completed"},
        {"component": "optional_abc_comparator", "status": "NOT_RUN_DISABLED", "reason": "Optional external comparator was not required and no ABC computation was performed"},
        {"component": "optional_differential_validation", "status": "NOT_AVAILABLE_INPUT", "reason": "DEFERRED_NO_VALID_COMPACT_QUANTITATIVE_INPUT"},
    ]
    write_tsv(out / "external_validation_status.tsv", ["component", "status", "reason"], statuses)

    metrics = read_tsv(Path(args.crispr) / "crispr_benchmark_metrics.tsv")
    distance_rows = [r for r in metrics if r["source_study"] == "ALL" and r["subset"] == "ALL_DISTAL_ELIGIBLE" and r["method"].startswith("candidate_")]
    distance_plot = [{"distance_kb": str(int(r["max_distance"]) // 1000), "precision": r["precision"], "recall": r["recall"], "F1": r["F1"]} for r in distance_rows]
    write_tsv(out / "figures" / "crispr_distance.tsv", ["distance_kb", "precision", "recall", "F1"], distance_plot)
    bar_svg(out / "figures" / "crispr_distance.svg", "CRISPR benchmark F1 by frozen distance", distance_plot, "distance_kb", "F1", 1)

    ccre = read_tsv(Path(args.ccre) / "ccre_concordance_summary.tsv")
    ccre_any = [{"evidence_class": r["evidence_class"], "overlap_proportion": r["overlap_proportion"], "total_regions": r["total_regions"]} for r in ccre if r["ccre_category"] == "ANY_CCRE"]
    write_tsv(out / "figures" / "ccre_overlap.tsv", ["evidence_class", "overlap_proportion", "total_regions"], ccre_any)
    bar_svg(out / "figures" / "ccre_overlap.svg", "ENCODE K562 cCRE annotation concordance", ccre_any, "evidence_class", "overlap_proportion", 1)

    plausibility = read_tsv(Path(args.plausibility) / "biological_plausibility_summary.tsv")
    plot_rows = [{"expression_stratum": r["expression_stratum"], "promoter_atac_proportion": r["promoter_atac_proportion"], "promoter_h3k27ac_proportion": r["promoter_h3k27ac_proportion"]} for r in plausibility]
    write_tsv(out / "figures" / "plausibility_promoter_evidence.tsv", ["expression_stratum", "promoter_atac_proportion", "promoter_h3k27ac_proportion"], plot_rows)
    bar_svg(out / "figures" / "plausibility_promoter_evidence.svg", "Promoter ATAC evidence by expression stratum", plot_rows, "expression_stratum", "promoter_atac_proportion", 1)

    important = [Path(args.sources), Path(args.candidate_freeze),
                 Path(args.ccre) / "ccre_concordance_summary.tsv",
                 Path(args.crispr) / "crispr_benchmark_pairs.tsv",
                 Path(args.crispr) / "crispr_benchmark_metrics.tsv",
                 Path(args.crispr) / "crispr_failure_classes.tsv",
                 Path(args.plausibility) / "biological_plausibility_summary.tsv"]
    checksums = {path.name: sha256(path) for path in important}
    provenance = {
        "schema_name": "external_validation_provenance", "schema_version": "1.0.0",
        "project": "regulatory-multiomics-pipeline", "project_version": "0.1.0",
        "author": "Sarvesh Galgale", "organism": "Homo sapiens", "cell_type": "K562",
        "genome_build": "GRCh38 primary assembly", "annotation": "GENCODE release 50",
        "coordinate_system": "0-based-half-open", "allowed_chromosomes": [str(i) for i in range(1, 23)] + ["X", "Y"],
        "integration_match_type": "REFERENCE_COMPARABLE",
        "primary_distal_distance_bp": 100000, "sensitivity_distances_bp": [50000, 100000, 250000],
        "promoter_definition": {"upstream_bp": 2000, "downstream_bp": 500, "strand_aware": True},
        "external_element_overlap": "strict positive base-pair intersection",
        "midpoint_mapping_sensitivity": "external element midpoint contained in validation region",
        "positive_label": "Regulated=TRUE: adjusted p-value below 0.05 and negative expression effect",
        "negative_label": "Regulated=FALSE among source-defined valid connections",
        "prediction_frozen_before_truth_comparison": True,
        "ccre_independence_class": "ANNOTATION_CONCORDANCE_NONINDEPENDENT",
        "crispr_independence_class": "FUNCTIONAL_INDEPENDENT",
        "optional_abc_comparator": "NOT_RUN_DISABLED",
        "optional_differential_validation": "DEFERRED_NO_VALID_COMPACT_QUANTITATIVE_INPUT",
        "checksums": checksums, "component_statuses": statuses,
        "software": {"Python": ".".join(map(str, sys.version_info[:3])), "interval_engine": "deterministic Python standard library"},
    }
    (out / "external_validation_provenance.json").write_text(json.dumps(provenance, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    print(f"statuses={len(statuses)} figures=3 provenance_inputs={len(checksums)}")


if __name__ == "__main__":
    main()
