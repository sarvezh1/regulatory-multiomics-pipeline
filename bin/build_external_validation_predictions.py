#!/usr/bin/env python3
"""Build frozen-rule full-genome regulatory candidates for external validation."""

import argparse
import bisect
import csv
import hashlib
from collections import defaultdict
from pathlib import Path


ALLOWED = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
ORDER = {chromosome: index for index, chromosome in enumerate(ALLOWED)}


def read_tsv(path):
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, fields, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def interval_index(rows, start_field="start", end_field="end"):
    grouped = defaultdict(list)
    maximum = 0
    for row in rows:
        start, end = int(row[start_field]), int(row[end_field])
        grouped[row["chromosome"]].append((start, end, row))
        maximum = max(maximum, end - start)
    result = {}
    for chromosome, items in grouped.items():
        items.sort(key=lambda item: (item[0], item[1]))
        result[chromosome] = ([item[0] for item in items], items)
    return result, maximum


def hits(index, maximum_length, chromosome, start, end):
    if chromosome not in index:
        return []
    starts, items = index[chromosome]
    low = bisect.bisect_left(starts, start - maximum_length)
    high = bisect.bisect_left(starts, end)
    return [row for item_start, item_end, row in items[low:high] if item_end > start]


def overlap(start_a, end_a, start_b, end_b):
    return max(0, min(end_a, end_b) - max(start_a, start_b))


def sorted_rows(rows):
    return sorted(
        rows,
        key=lambda row: (
            ORDER[row["chromosome"]],
            int(row.get("start", row.get("overlap_start", 0))),
            int(row.get("end", row.get("overlap_end", 0))),
            row.get("region_id", ""),
        ),
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tss", required=True)
    parser.add_argument("--promoters", required=True)
    parser.add_argument("--atac", required=True)
    parser.add_argument("--h3k27ac", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--distances", nargs="+", type=int, default=[50000, 100000, 250000])
    args = parser.parse_args()

    output = Path(args.outdir)
    tss = read_tsv(args.tss)
    promoters = read_tsv(args.promoters)
    atac = read_tsv(args.atac)
    chip = read_tsv(args.h3k27ac)
    for rows, start_field, end_field in (
        (promoters, "promoter_start", "promoter_end"),
        (atac, "start", "end"),
        (chip, "start", "end"),
    ):
        for row in rows:
            row["start"] = int(row[start_field])
            row["end"] = int(row[end_field])

    promoter_index, promoter_max = interval_index(promoters)
    chip_index, chip_max = interval_index(chip)
    atac_index, atac_max = interval_index(atac)

    promoter_atac = []
    promoter_chip = []
    promoter_atac_ids = set()
    promoter_chip_ids = set()
    for region, destination, observed in (
        (atac, promoter_atac, promoter_atac_ids),
        (chip, promoter_chip, promoter_chip_ids),
    ):
        for source in region:
            for promoter in hits(
                promoter_index,
                promoter_max,
                source["chromosome"],
                source["start"],
                source["end"],
            ):
                overlap_bp = overlap(
                    source["start"], source["end"], promoter["start"], promoter["end"]
                )
                if overlap_bp:
                    observed.add(source["region_id"])
                    destination.append(
                        {
                            "source_region_id": source["region_id"],
                            "promoter_id": promoter["promoter_id"],
                            "gene_id_base": promoter["gene_id_base"],
                            "chromosome": source["chromosome"],
                            "overlap_bp": overlap_bp,
                        }
                    )

    colocated = []
    atac_to_chip = defaultdict(set)
    chip_to_atac = defaultdict(set)
    for atac_region in atac:
        for chip_region in hits(
            chip_index,
            chip_max,
            atac_region["chromosome"],
            atac_region["start"],
            atac_region["end"],
        ):
            start = max(atac_region["start"], chip_region["start"])
            end = min(atac_region["end"], chip_region["end"])
            if end <= start:
                continue
            token = hashlib.sha256(
                f"{atac_region['region_id']}|{chip_region['region_id']}|{start}|{end}".encode()
            ).hexdigest()[:16]
            promoter_overlap = bool(
                hits(promoter_index, promoter_max, atac_region["chromosome"], start, end)
            )
            colocated.append(
                {
                    "colocalized_region_id": f"COLOC_{token}",
                    "atac_region_id": atac_region["region_id"],
                    "h3k27ac_region_id": chip_region["region_id"],
                    "chromosome": atac_region["chromosome"],
                    "overlap_start": start,
                    "overlap_end": end,
                    "overlap_bp": end - start,
                    "promoter_overlap": str(promoter_overlap).lower(),
                    "relationship": "ATAC_H3K27AC_COLOCATED",
                }
            )
            atac_to_chip[atac_region["region_id"]].add(chip_region["region_id"])
            chip_to_atac[chip_region["region_id"]].add(atac_region["region_id"])

    tss_by_chromosome = defaultdict(list)
    for row in tss:
        if row["chromosome"] in ORDER:
            tss_by_chromosome[row["chromosome"]].append((int(row["tss_start"]), row))
    for chromosome in tss_by_chromosome:
        tss_by_chromosome[chromosome].sort(key=lambda item: (item[0], item[1]["tss_id"]))
    tss_positions = {
        chromosome: [entry[0] for entry in entries]
        for chromosome, entries in tss_by_chromosome.items()
    }
    tss_rows_by_position = {}
    for chromosome, entries in tss_by_chromosome.items():
        grouped = defaultdict(list)
        for position, row in entries:
            grouped[position].append(row)
        tss_rows_by_position[chromosome] = grouped

    distal_atac = [row for row in atac if row["region_id"] not in promoter_atac_ids]
    prediction_sets = {}
    for maximum_distance in sorted(set(args.distances)):
        predictions = []
        for region in distal_atac:
            entries = tss_by_chromosome.get(region["chromosome"], [])
            if not entries:
                continue
            positions = tss_positions[region["chromosome"]]
            insertion = bisect.bisect_left(positions, region["start"])
            candidate_positions = []
            if insertion:
                candidate_positions.append(positions[insertion - 1])
            right = bisect.bisect_left(positions, region["end"])
            if right < len(positions):
                candidate_positions.append(positions[right])
            if not candidate_positions:
                continue
            distances = {
                position: (
                    region["start"] - position
                    if position < region["start"]
                    else position - (region["end"] - 1)
                )
                for position in candidate_positions
            }
            nearest_distance = min(distances.values())
            if nearest_distance > maximum_distance:
                continue
            nearest_positions = {
                position for position, distance in distances.items() if distance == nearest_distance
            }
            for position in sorted(nearest_positions):
                for tss_row in tss_rows_by_position[region["chromosome"]][position]:
                    predictions.append(
                        {
                            "region_id": region["region_id"],
                            "candidate_gene_id": tss_row["gene_id_base"],
                            "candidate_tss_id": tss_row["tss_id"],
                            "distance_to_tss": nearest_distance,
                            "link_method": "NEAREST_TSS_CANDIDATE",
                            "max_distance": maximum_distance,
                            "source_modality": "ATAC",
                            "anchor_type": "ATAC_ANCHORED",
                            "h3k27ac_region_ids": ";".join(
                                sorted(atac_to_chip.get(region["region_id"], set()))
                            )
                            or "NA",
                        }
                    )
        predictions.sort(
            key=lambda row: (
                row["region_id"], row["candidate_gene_id"], row["candidate_tss_id"]
            )
        )
        prediction_sets[maximum_distance] = predictions
        write_tsv(
            output / "predictions" / f"candidate_links.{maximum_distance}.tsv",
            [
                "region_id",
                "candidate_gene_id",
                "candidate_tss_id",
                "distance_to_tss",
                "link_method",
                "max_distance",
                "source_modality",
                "anchor_type",
                "h3k27ac_region_ids",
            ],
            predictions,
        )

    validation_regions = []
    for region in atac:
        if region["region_id"] in promoter_atac_ids:
            evidence_class = "PROMOTER_ATAC"
        elif region["region_id"] in atac_to_chip:
            evidence_class = "DISTAL_ATAC_H3K27AC_COLOCATED"
        else:
            evidence_class = "DISTAL_ATAC_ONLY"
        validation_regions.append({**region, "evidence_class": evidence_class})
    for region in chip:
        if region["region_id"] in promoter_chip_ids:
            evidence_class = "PROMOTER_H3K27AC"
        elif region["region_id"] not in chip_to_atac:
            evidence_class = "DISTAL_H3K27AC_ONLY"
        else:
            continue
        validation_regions.append({**region, "evidence_class": evidence_class})

    write_tsv(
        output / "relationships" / "promoter_atac.tsv",
        ["source_region_id", "promoter_id", "gene_id_base", "chromosome", "overlap_bp"],
        promoter_atac,
    )
    write_tsv(
        output / "relationships" / "promoter_h3k27ac.tsv",
        ["source_region_id", "promoter_id", "gene_id_base", "chromosome", "overlap_bp"],
        promoter_chip,
    )
    write_tsv(
        output / "relationships" / "atac_h3k27ac_colocalization.tsv",
        [
            "colocalized_region_id",
            "atac_region_id",
            "h3k27ac_region_id",
            "chromosome",
            "overlap_start",
            "overlap_end",
            "overlap_bp",
            "promoter_overlap",
            "relationship",
        ],
        colocated,
    )
    write_tsv(
        output / "regions" / "validation_regions.tsv",
        ["region_id", "chromosome", "start", "end", "coordinate_system", "evidence_class"],
        sorted_rows(validation_regions),
    )
    print(
        " ".join(
            [
                f"promoter_atac_links={len(promoter_atac)}",
                f"promoter_h3k27ac_links={len(promoter_chip)}",
                f"colocalizations={len(colocated)}",
                *[
                    f"candidates_{distance}={len(prediction_sets[distance])}"
                    for distance in sorted(prediction_sets)
                ],
            ]
        )
    )


if __name__ == "__main__":
    main()
