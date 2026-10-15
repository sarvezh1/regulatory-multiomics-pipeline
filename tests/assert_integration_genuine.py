#!/usr/bin/env python3
import csv
import json
import pathlib


ROOT = pathlib.Path("tests/output/integration_k562/integration")


def rows(relative):
    with (ROOT / relative).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


required = [
    "compatibility/integration.contract_compatibility.tsv",
    "evidence/promoter_evidence.tsv",
    "evidence/promoter_methylation_by_sample.tsv",
    "relationships/atac_promoter.tsv",
    "relationships/h3k27ac_promoter.tsv",
    "relationships/cpg_promoter.tsv",
    "relationships/atac_h3k27ac_colocalization.tsv",
    "relationships/distal_candidate_links.tsv",
    "provenance/integration.analysis_status.tsv",
    "provenance/integration.provenance.json",
]
for relative in required:
    path = ROOT / relative
    assert path.is_file() and path.stat().st_size > 0, relative

compatibility = rows("compatibility/integration.contract_compatibility.tsv")
assert {row["modality"] for row in compatibility} == {"RNA", "ATAC", "WGBS", "H3K27ac"}
assert all(row["reference_compatibility"] == "REFERENCE_COMPATIBLE" for row in compatibility)

statuses = {row["component"]: row["status"] for row in rows("provenance/integration.analysis_status.tsv")}
assert statuses == {
    "integration_reference": "COMPLETE",
    "integration_promoter": "COMPLETE",
    "integration_differential": "NOT_RUN_DISABLED",
    "integration_distal": "COMPLETE",
}

promoters = rows("evidence/promoter_evidence.tsv")
assert promoters and all(row["integration_match_type"] == "REFERENCE_COMPARABLE" for row in promoters)
for row in promoters:
    assert row["tested_modality_count"] == "0"
    assert row["rna_state"] == "NOT_TESTED"
    assert row["atac_state"] == "NOT_TESTED"
    assert row["methylation_state"] == "NOT_TESTED"
    assert row["h3k27ac_state"] == "NOT_TESTED"

methylation = rows("evidence/promoter_methylation_by_sample.tsv")
assert methylation
assert all(row["evidence_type"] == "DESCRIPTIVE_PERCENTAGE_WITH_COVERAGE" for row in methylation)
assert all(row["total_methylated_reads"] == "NA" for row in methylation)
assert all(row["total_unmethylated_reads"] == "NA" for row in methylation)
assert all(0.0 <= float(row["coverage_weighted_methylation_fraction"]) <= 1.0 for row in methylation)

cpg_links = rows("relationships/cpg_promoter.tsv")
assert cpg_links
assert len({(row["cpg_id"], row["promoter_id"], row["sample_id"]) for row in cpg_links}) == len(cpg_links)
assert all(int(row["overlap_bp"]) > 0 for row in cpg_links)

assert rows("relationships/atac_promoter.tsv")
assert rows("relationships/h3k27ac_promoter.tsv")
assert rows("relationships/atac_h3k27ac_colocalization.tsv")
assert rows("relationships/distal_candidate_links.tsv")

with (ROOT / "provenance/integration.provenance.json").open(encoding="utf-8") as handle:
    provenance = json.load(handle)
assert provenance["integration_match_type"] == "REFERENCE_COMPARABLE"
assert provenance["contrast_compatibility"] == "CONTRAST_UNAVAILABLE"
assert provenance["no_cross_assay_pvalue_combination"] is True
assert provenance["no_combined_regulatory_score"] is True

print(
    "K562 descriptive integration assertions passed: "
    f"{len(promoters)} promoters, {len(cpg_links)} CpG-promoter links"
)
