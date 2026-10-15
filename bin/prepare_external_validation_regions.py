#!/usr/bin/env python3
"""Prepare full-genome K562 interval contracts for external validation."""

import argparse
import csv
import gzip
import hashlib
from pathlib import Path


ALLOWED = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}


def write_tsv(path, fieldnames, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def convert_peaks(path, accession, prefix):
    unique = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            columns = line.rstrip("\n").split("\t")
            if len(columns) < 3:
                raise ValueError(f"{path.name}: malformed BED row {line_number}")
            chromosome = columns[0]
            if chromosome not in ALLOWED:
                continue
            start, end = int(columns[1]), int(columns[2])
            if start < 0 or end <= start:
                raise ValueError(f"{path.name}: invalid interval at row {line_number}")
            token = hashlib.sha256(
                f"{chromosome}:{start}-{end}".encode("utf-8")
            ).hexdigest()[:16]
            unique[(chromosome, start, end)] = {
                "region_id": f"{prefix}_{accession}_{token}",
                "chromosome": chromosome,
                "start": start,
                "end": end,
                "coordinate_system": "0-based-half-open",
            }
    chromosome_order = {f"chr{i}": i for i in range(1, 23)} | {
        "chrX": 23,
        "chrY": 24,
    }
    return sorted(
        unique.values(),
        key=lambda row: (
            chromosome_order[row["chromosome"]],
            row["start"],
            row["end"],
            row["region_id"],
        ),
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--atac", required=True)
    parser.add_argument("--h3k27ac", required=True)
    parser.add_argument("--rna", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    output = Path(args.output_dir)
    atac = convert_peaks(Path(args.atac), "ENCFF057UYP", "ATAC")
    chip = convert_peaks(Path(args.h3k27ac), "ENCFF038DDS", "H3K27AC")
    fields = ["region_id", "chromosome", "start", "end", "coordinate_system"]
    write_tsv(output / "k562.atac_regions.full_genome.tsv", fields, atac)
    write_tsv(output / "k562.h3k27ac_regions.full_genome.tsv", fields, chip)

    manifest_fields = [
        "modality",
        "data_type",
        "path",
        "source_genome_build",
        "coordinate_normalization_status",
        "source_annotation",
        "status",
        "contrast_variable",
        "numerator",
        "denominator",
        "integration_match_type",
    ]
    common = {
        "source_genome_build": "GRCh38",
        "coordinate_normalization_status": "NORMALIZED",
        "source_annotation": "NOT_APPLICABLE",
        "status": "COMPLETE",
        "contrast_variable": "NA",
        "numerator": "NA",
        "denominator": "NA",
        "integration_match_type": "REFERENCE_COMPARABLE",
    }
    rows = [
        {
            **common,
            "modality": "RNA",
            "data_type": "rna_baseline",
            "path": Path(args.rna).as_posix(),
            "coordinate_normalization_status": "NOT_APPLICABLE",
            "source_annotation": "GENCODE V29 source quantification; identifiers retained",
        },
        {
            **common,
            "modality": "ATAC",
            "data_type": "atac_regions",
            "path": (output / "k562.atac_regions.full_genome.tsv").as_posix(),
        },
        {
            **common,
            "modality": "H3K27ac",
            "data_type": "h3k27ac_regions",
            "path": (output / "k562.h3k27ac_regions.full_genome.tsv").as_posix(),
        },
    ]
    write_tsv(output / "manifest.k562_full_genome.tsv", manifest_fields, rows)
    print(f"atac_regions={len(atac)} h3k27ac_regions={len(chip)}")


if __name__ == "__main__":
    main()
