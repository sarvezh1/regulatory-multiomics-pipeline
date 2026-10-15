#!/usr/bin/env python3
"""Measure concordance between frozen K562 region classes and ENCODE cCREs."""

import argparse
import csv
import gzip
from collections import Counter, defaultdict
from pathlib import Path


ALLOWED = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]


def read_regions(path):
    grouped = defaultdict(list)
    with Path(path).open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            row["start"] = int(row["start"])
            row["end"] = int(row["end"])
            grouped[row["chromosome"]].append(row)
    for chromosome in grouped:
        grouped[chromosome].sort(key=lambda row: (row["start"], row["end"], row["region_id"]))
    return grouped


def write_tsv(path, fields, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--regions", required=True)
    parser.add_argument("--ccre", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    regions = read_regions(args.regions)
    total = Counter(
        row["evidence_class"] for chromosome in regions.values() for row in chromosome
    )
    overlaps = []
    region_any = defaultdict(set)
    region_category = defaultdict(set)
    category_vocabulary = Counter()
    ccre_rows = 0

    current_chromosome = None
    current_regions = []
    pointer = 0
    active = []
    with gzip.open(args.ccre, "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                raise ValueError(f"cCRE row {line_number} has fewer than 11 fields")
            chromosome = fields[0]
            if chromosome not in regions:
                continue
            start, end = int(fields[1]), int(fields[2])
            if start < 0 or end <= start:
                raise ValueError(f"cCRE row {line_number} has invalid coordinates")
            category = fields[9]
            category_vocabulary[category] += 1
            ccre_rows += 1
            if chromosome != current_chromosome:
                current_chromosome = chromosome
                current_regions = regions[chromosome]
                pointer = 0
                active = []
            while pointer < len(current_regions) and current_regions[pointer]["start"] < end:
                active.append(current_regions[pointer])
                pointer += 1
            active = [region for region in active if region["end"] > start]
            for region in active:
                overlap_bp = max(0, min(end, region["end"]) - max(start, region["start"]))
                if not overlap_bp:
                    continue
                evidence_class = region["evidence_class"]
                region_any[evidence_class].add(region["region_id"])
                region_category[(evidence_class, category)].add(region["region_id"])
                overlaps.append(
                    {
                        "region_id": region["region_id"],
                        "evidence_class": evidence_class,
                        "ccre_id": fields[3],
                        "ccre_category": category,
                        "ccre_classification_set": fields[10],
                        "chromosome": chromosome,
                        "overlap_bp": overlap_bp,
                    }
                )

    summary = []
    for evidence_class in sorted(total):
        count = total[evidence_class]
        any_count = len(region_any[evidence_class])
        summary.append(
            {
                "evidence_class": evidence_class,
                "ccre_category": "ANY_CCRE",
                "total_regions": count,
                "overlapping_regions": any_count,
                "overlap_proportion": any_count / count if count else "NA",
            }
        )
        for category in sorted(category_vocabulary):
            category_count = len(region_category[(evidence_class, category)])
            summary.append(
                {
                    "evidence_class": evidence_class,
                    "ccre_category": category,
                    "total_regions": count,
                    "overlapping_regions": category_count,
                    "overlap_proportion": category_count / count if count else "NA",
                }
            )

    write_tsv(
        Path(args.outdir) / "ccre_region_overlaps.tsv",
        [
            "region_id",
            "evidence_class",
            "ccre_id",
            "ccre_category",
            "ccre_classification_set",
            "chromosome",
            "overlap_bp",
        ],
        overlaps,
    )
    write_tsv(
        Path(args.outdir) / "ccre_concordance_summary.tsv",
        [
            "evidence_class",
            "ccre_category",
            "total_regions",
            "overlapping_regions",
            "overlap_proportion",
        ],
        summary,
    )
    write_tsv(
        Path(args.outdir) / "ccre_category_vocabulary.tsv",
        ["ccre_category", "source_record_count"],
        [
            {"ccre_category": category, "source_record_count": count}
            for category, count in sorted(category_vocabulary.items())
        ],
    )
    print(
        f"ccre_primary_records={ccre_rows} overlaps={len(overlaps)} "
        f"region_classes={len(total)} categories={len(category_vocabulary)}"
    )


if __name__ == "__main__":
    main()
