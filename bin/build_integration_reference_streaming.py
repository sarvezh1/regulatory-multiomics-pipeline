#!/usr/bin/env python3
"""Derive the frozen integration TSS/promoter reference with bounded memory."""

import argparse
import csv
import gzip
import re
from pathlib import Path


ALLOWED = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
ORDER = {chromosome: index for index, chromosome in enumerate(ALLOWED)}
ATTRIBUTE = re.compile(r'([A-Za-z0-9_]+) "([^"]*)"')


def opener(path):
    return gzip.open(path, "rt", encoding="utf-8") if str(path).endswith(".gz") else open(path, encoding="utf-8")


def fasta_lengths(path):
    lengths = {}
    current = None
    with opener(path) as handle:
        for line in handle:
            if line.startswith(">"):
                current = line[1:].split()[0]
                if current in ALLOWED:
                    lengths[current] = 0
                continue
            if current in lengths:
                lengths[current] += len(line.strip())
    if not lengths:
        raise ValueError("FASTA contains no allowed primary chromosome")
    return lengths


def write_tsv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("gtf")
    parser.add_argument("fasta")
    parser.add_argument("outdir")
    parser.add_argument("annotation")
    args = parser.parse_args()
    lengths = fasta_lengths(args.fasta)
    groups = {}
    with opener(args.gtf) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line or line.startswith("#"):
                continue
            columns = line.rstrip("\n").split("\t")
            if len(columns) != 9:
                raise ValueError(f"GTF row {line_number} does not have 9 columns")
            chromosome, _, feature, start, end, _, strand, _, attributes = columns
            if feature != "transcript" or chromosome not in ORDER:
                continue
            parsed = dict(ATTRIBUTE.findall(attributes))
            if "gene_id" not in parsed or "transcript_id" not in parsed:
                raise ValueError(f"GTF transcript row {line_number} lacks gene_id/transcript_id")
            gene_original = parsed["gene_id"]
            gene_base = re.sub(r"\.[0-9]+$", "", gene_original)
            position = int(start) if strand == "+" else int(end)
            key = (chromosome, position, strand, gene_base)
            group = groups.setdefault(
                key,
                {"gene_original": set(), "gene_name": set(), "transcript_id": set()},
            )
            group["gene_original"].add(gene_original)
            if parsed.get("gene_name"):
                group["gene_name"].add(parsed["gene_name"])
            group["transcript_id"].add(parsed["transcript_id"])
    missing = {key[0] for key in groups} - set(lengths)
    if missing:
        raise ValueError(f"GTF chromosomes absent from FASTA: {','.join(sorted(missing))}")

    tss = []
    for (chromosome, position, strand, gene_base), values in groups.items():
        tss.append(
            {
                "tss_id": f"TSS_{gene_base}_{chromosome}_{position}_{'PLUS' if strand == '+' else 'MINUS'}",
                "gene_id_original": ";".join(sorted(values["gene_original"])),
                "gene_id_base": gene_base,
                "gene_name": ";".join(sorted(values["gene_name"])) or "NA",
                "transcript_ids": ";".join(sorted(values["transcript_id"])),
                "chromosome": chromosome,
                "tss_start": position - 1,
                "tss_end": position,
                "strand": strand,
                "coordinate_system": "0-based-half-open",
                "source_annotation": args.annotation,
            }
        )
    tss.sort(key=lambda row: (ORDER[row["chromosome"]], row["tss_start"], row["strand"], row["gene_id_base"]))
    outdir = Path(args.outdir)
    tss_fields = list(tss[0])
    write_tsv(outdir / "integration.tss.tsv", tss_fields, tss)

    for upstream, downstream, label in (
        (2000, 500, "m2000_p500"),
        (1000, 500, "m1000_p500"),
        (2000, 1000, "m2000_p1000"),
    ):
        promoters = []
        for row in tss:
            position = int(row["tss_end"])
            if row["strand"] == "+":
                start = max(0, position - upstream - 1)
                end = min(lengths[row["chromosome"]], position + downstream - 1)
            else:
                start = max(0, position - downstream)
                end = min(lengths[row["chromosome"]], position + upstream)
            promoters.append(
                {
                    "promoter_id": f"PROM_{row['tss_id']}_U{upstream}_D{downstream}",
                    "tss_id": row["tss_id"],
                    "gene_id_original": row["gene_id_original"],
                    "gene_id_base": row["gene_id_base"],
                    "gene_name": row["gene_name"],
                    "chromosome": row["chromosome"],
                    "promoter_start": start,
                    "promoter_end": end,
                    "strand": row["strand"],
                    "coordinate_system": "0-based-half-open",
                    "upstream_bp": upstream,
                    "downstream_bp": downstream,
                    "source_annotation": args.annotation,
                }
            )
        write_tsv(outdir / f"integration.promoters.{label}.tsv", list(promoters[0]), promoters)
    (outdir / "integration_reference.versions.txt").write_text(
        "runtime\tPython standard library\nalgorithm\tintegration all-distinct-transcript-TSS policy\n",
        encoding="utf-8",
    )
    print(f"tss={len(tss)} allowed_chromosomes={len(lengths)}")


if __name__ == "__main__":
    main()
