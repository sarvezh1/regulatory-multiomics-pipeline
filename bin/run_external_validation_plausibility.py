#!/usr/bin/env python3
"""Run descriptive K562 cross-gene plausibility summaries."""

import argparse
import csv
import math
import statistics
from collections import defaultdict
from pathlib import Path


NA_VALUES = {"", "NA", "NaN", "nan"}


def read_tsv(path):
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, fields, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def number(value):
    if value in NA_VALUES or value is None:
        return None
    return float(value)


def quantile(values, probability):
    ordered = sorted(values)
    if not ordered:
        return None
    position = (len(ordered) - 1) * probability
    low = math.floor(position)
    high = math.ceil(position)
    if low == high:
        return ordered[low]
    fraction = position - low
    return ordered[low] * (1 - fraction) + ordered[high] * fraction


def present(value):
    return "NA" if value is None else value


def ranks(values):
    order = sorted(range(len(values)), key=lambda index: values[index])
    result = [0.0] * len(values)
    cursor = 0
    while cursor < len(order):
        end = cursor + 1
        while end < len(order) and values[order[end]] == values[order[cursor]]:
            end += 1
        rank = (cursor + 1 + end) / 2
        for index in order[cursor:end]:
            result[index] = rank
        cursor = end
    return result


def pearson(x, y):
    if len(x) < 2:
        return None
    mean_x, mean_y = statistics.mean(x), statistics.mean(y)
    numerator = sum((a - mean_x) * (b - mean_y) for a, b in zip(x, y))
    denominator = math.sqrt(
        sum((a - mean_x) ** 2 for a in x) * sum((b - mean_y) ** 2 for b in y)
    )
    return numerator / denominator if denominator else None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--promoter-evidence", required=True)
    parser.add_argument("--promoter-methylation", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    promoters = read_tsv(args.promoter_evidence)
    methylation_rows = read_tsv(args.promoter_methylation)
    methylation = defaultdict(list)
    for row in methylation_rows:
        fraction = number(row.get("coverage_weighted_methylation_fraction"))
        coverage = number(row.get("total_coverage"))
        if fraction is not None and coverage is not None and coverage > 0:
            methylation[row["promoter_id"]].append((fraction, coverage))

    by_gene = defaultdict(list)
    for row in promoters:
        expression = number(row.get("rna_descriptive_expression"))
        if expression is not None:
            by_gene[row["gene_id_base"]].append(row)
    expression_values = []
    gene_expression = {}
    for gene, rows in by_gene.items():
        values = sorted(
            {
                number(row["rna_descriptive_expression"])
                for row in rows
                if number(row["rna_descriptive_expression"]) is not None
            }
        )
        if len(values) != 1:
            raise ValueError(f"gene {gene} has inconsistent descriptive RNA values")
        gene_expression[gene] = values[0]
        expression_values.append(values[0])
    q25 = quantile(expression_values, 0.25)
    q75 = quantile(expression_values, 0.75)

    gene_rows = []
    promoter_rows_by_stratum = defaultdict(list)
    for gene in sorted(by_gene):
        rows = by_gene[gene]
        expression = gene_expression[gene]
        if expression <= q25:
            stratum = "BOTTOM_QUARTILE"
        elif expression >= q75:
            stratum = "TOP_QUARTILE"
        else:
            stratum = "MIDDLE_50_PERCENT"
        atac_count = sum(int(row["n_atac_regions"]) > 0 for row in rows)
        chip_count = sum(int(row["n_h3k27ac_regions"]) > 0 for row in rows)
        weighted_sum = 0.0
        total_coverage = 0.0
        for row in rows:
            for fraction, coverage in methylation.get(row["promoter_id"], []):
                weighted_sum += fraction * coverage
                total_coverage += coverage
            promoter_rows_by_stratum[stratum].append(row)
        gene_rows.append(
            {
                "gene_id_base": gene,
                "gene_name": rows[0]["gene_name"],
                "rna_descriptive_expression": expression,
                "expression_stratum": stratum,
                "number_of_promoters": len(rows),
                "promoter_atac_fraction": atac_count / len(rows),
                "promoter_h3k27ac_fraction": chip_count / len(rows),
                "coverage_weighted_promoter_methylation": (
                    weighted_sum / total_coverage if total_coverage else "NA"
                ),
                "methylation_total_coverage": total_coverage if total_coverage else "NA",
                "integration_match_type": "REFERENCE_COMPARABLE",
            }
        )

    summary = []
    for stratum in ("BOTTOM_QUARTILE", "MIDDLE_50_PERCENT", "TOP_QUARTILE"):
        genes = [row for row in gene_rows if row["expression_stratum"] == stratum]
        promoter_rows = promoter_rows_by_stratum[stratum]
        methylation_values = []
        for row in promoter_rows:
            for fraction, _ in methylation.get(row["promoter_id"], []):
                methylation_values.append(fraction)
        summary.append(
            {
                "expression_stratum": stratum,
                "gene_count": len(genes),
                "promoter_count": len(promoter_rows),
                "promoters_with_atac": sum(
                    int(row["n_atac_regions"]) > 0 for row in promoter_rows
                ),
                "promoter_atac_proportion": (
                    sum(int(row["n_atac_regions"]) > 0 for row in promoter_rows)
                    / len(promoter_rows)
                    if promoter_rows
                    else "NA"
                ),
                "promoters_with_h3k27ac": sum(
                    int(row["n_h3k27ac_regions"]) > 0 for row in promoter_rows
                ),
                "promoter_h3k27ac_proportion": (
                    sum(int(row["n_h3k27ac_regions"]) > 0 for row in promoter_rows)
                    / len(promoter_rows)
                    if promoter_rows
                    else "NA"
                ),
                "promoters_with_methylation": len(methylation_values),
                "promoter_methylation_q25": present(quantile(methylation_values, 0.25)),
                "promoter_methylation_median": present(quantile(methylation_values, 0.5)),
                "promoter_methylation_q75": present(quantile(methylation_values, 0.75)),
                "integration_match_type": "REFERENCE_COMPARABLE",
            }
        )

    correlations = []
    for column in (
        "promoter_atac_fraction",
        "promoter_h3k27ac_fraction",
        "coverage_weighted_promoter_methylation",
    ):
        pairs = [
            (row["rna_descriptive_expression"], number(row[column]))
            for row in gene_rows
            if number(row[column]) is not None
        ]
        rho = pearson(ranks([x for x, _ in pairs]), ranks([y for _, y in pairs]))
        correlations.append(
            {
                "x_variable": "rna_descriptive_expression",
                "y_variable": column,
                "method": "SPEARMAN_DESCRIPTIVE",
                "rho": rho if rho is not None else "NA",
                "n_complete": len(pairs),
                "pvalue": "NOT_CALCULATED",
                "integration_match_type": "REFERENCE_COMPARABLE",
            }
        )

    output = Path(args.outdir)
    write_tsv(output / "biological_plausibility_gene_table.tsv", list(gene_rows[0]), gene_rows)
    write_tsv(output / "biological_plausibility_summary.tsv", list(summary[0]), summary)
    write_tsv(output / "biological_plausibility_correlations.tsv", list(correlations[0]), correlations)
    print(
        f"genes={len(gene_rows)} q25={q25:.6g} q75={q75:.6g} "
        "integration_match_type=REFERENCE_COMPARABLE"
    )


if __name__ == "__main__":
    main()
