#!/usr/bin/env python3
"""Benchmark frozen distal candidate links against tested K562 CRISPR pairs."""

import argparse
import bisect
import csv
import gzip
import math
from collections import Counter, defaultdict
from pathlib import Path


ALLOWED = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}
DISTANCE_BINS = ((0, 10_000, "0-10kb"), (10_000, 50_000, ">10-50kb"),
                 (50_000, 100_000, ">50-100kb"), (100_000, 250_000, ">100-250kb"),
                 (250_000, math.inf, ">250kb"))


def open_text(path):
    return gzip.open(path, "rt", encoding="utf-8-sig", newline="") if str(path).endswith(".gz") else open(path, encoding="utf-8-sig", newline="")


def read_tsv(path):
    with open_text(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, fields, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def interval_index(rows, chrom="chromosome", start="start", end="end"):
    grouped = defaultdict(list)
    maximum = 0
    for row in rows:
        a, b = int(row[start]), int(row[end])
        grouped[row[chrom]].append((a, b, row))
        maximum = max(maximum, b - a)
    result = {}
    for chromosome, values in grouped.items():
        values.sort(key=lambda item: (item[0], item[1]))
        result[chromosome] = ([item[0] for item in values], values)
    return result, maximum


def overlaps(index, maximum, chromosome, start, end):
    if chromosome not in index:
        return []
    starts, rows = index[chromosome]
    lo = bisect.bisect_left(starts, start - maximum)
    hi = bisect.bisect_left(starts, end)
    return [row for a, b, row in rows[lo:hi] if b > start]


def safe_div(a, b):
    return a / b if b else math.nan


def fmt(value):
    if isinstance(value, float):
        return "NA" if math.isnan(value) else f"{value:.12g}"
    return str(value)


def logchoose(n, k):
    if k < 0 or k > n:
        return -math.inf
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def fisher(tp, fp, fn, tn):
    row1, row2, col1, total = tp + fp, fn + tn, tp + fn, tp + fp + fn + tn
    if total == 0:
        return math.nan, math.nan, math.nan, math.nan
    lo, hi = max(0, col1 - row2), min(row1, col1)
    def probability(x):
        return math.exp(logchoose(row1, x) + logchoose(row2, col1 - x) - logchoose(total, col1))
    observed = probability(tp)
    pvalue = min(1.0, sum(probability(x) for x in range(lo, hi + 1) if probability(x) <= observed * (1 + 1e-12)))
    odds = safe_div(tp * tn, fp * fn)
    if min(tp, fp, fn, tn) > 0:
        se = math.sqrt(1 / tp + 1 / fp + 1 / fn + 1 / tn)
        low = math.exp(math.log(odds) - 1.96 * se)
        high = math.exp(math.log(odds) + 1.96 * se)
    else:
        low = high = math.nan
    return odds, pvalue, low, high


def metric_row(method, distance, subset, records, source_study="ALL"):
    tp = sum(r["external_positive"] and r[method] for r in records)
    fp = sum((not r["external_positive"]) and r[method] for r in records)
    tn = sum((not r["external_positive"]) and (not r[method]) for r in records)
    fn = sum(r["external_positive"] and (not r[method]) for r in records)
    precision = safe_div(tp, tp + fp)
    recall = safe_div(tp, tp + fn)
    specificity = safe_div(tn, tn + fp)
    f1 = safe_div(2 * precision * recall, precision + recall) if not (math.isnan(precision) or math.isnan(recall)) else math.nan
    balanced = (recall + specificity) / 2 if not (math.isnan(recall) or math.isnan(specificity)) else math.nan
    mcc_denom = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = safe_div(tp * tn - fp * fn, mcc_denom)
    odds, fisher_p, odds_low, odds_high = fisher(tp, fp, fn, tn)
    return {k: fmt(v) for k, v in {
        "method": method, "max_distance": distance, "subset": subset, "source_study": source_study,
        "n_scored": len(records), "positive_records": tp + fn, "negative_records": tn + fp,
        "TP": tp, "FP": fp, "TN": tn, "FN": fn, "precision": precision, "recall": recall,
        "specificity": specificity, "F1": f1, "balanced_accuracy": balanced, "MCC": mcc,
        "fisher_odds_ratio": odds, "fisher_two_sided_pvalue": fisher_p,
        "odds_ratio_ci95_low_approx": odds_low, "odds_ratio_ci95_high_approx": odds_high,
    }.items()}


def nearest_tss(tss_by_chr, chromosome, start, end):
    entries = tss_by_chr.get(chromosome, [])
    if not entries:
        return math.inf, set()
    positions = [entry[0] for entry in entries]
    left = bisect.bisect_left(positions, start)
    candidates = []
    if left:
        candidates.append(positions[left - 1])
    right = bisect.bisect_left(positions, end)
    if right < len(positions):
        candidates.append(positions[right])
    distances = {p: start - p if p < start else p - (end - 1) for p in candidates}
    best = min(distances.values())
    genes = {row["gene_id_base"] for p, row in entries if p in distances and distances[p] == best}
    return best, genes


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmark", required=True)
    parser.add_argument("--tss", required=True)
    parser.add_argument("--promoters", required=True)
    parser.add_argument("--regions", required=True)
    parser.add_argument("--colocalization", required=True)
    parser.add_argument("--prediction", action="append", nargs=2, metavar=("DISTANCE", "PATH"), required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    out = Path(args.outdir)
    source = read_tsv(args.benchmark)
    tss = read_tsv(args.tss)
    promoters = read_tsv(args.promoters)
    regions = read_tsv(args.regions)
    coloc = read_tsv(args.colocalization)
    distances = sorted(int(distance) for distance, _ in args.prediction)

    genes = defaultdict(set)
    gene_symbols = defaultdict(set)
    tss_by_chr = defaultdict(list)
    for row in tss:
        genes[row["gene_id_base"]].add(row["gene_name"])
        gene_symbols[row["gene_name"]].add(row["gene_id_base"])
        tss_by_chr[row["chromosome"]].append((int(row["tss_start"]), row))
    for chromosome in tss_by_chr:
        tss_by_chr[chromosome].sort(key=lambda item: (item[0], item[1]["tss_id"]))

    promoters_by_gene = defaultdict(list)
    for row in promoters:
        promoters_by_gene[row["gene_id_base"]].append((int(row["promoter_start"]), int(row["promoter_end"]), row["chromosome"]))

    region_index, region_max = interval_index(regions)
    colocated_atac = {row["atac_region_id"] for row in coloc if row["promoter_overlap"].lower() == "false"}
    predictions = {}
    prediction_regions = {}
    for distance, path in args.prediction:
        entries = read_tsv(path)
        key = int(distance)
        predictions[key] = {(r["region_id"], r["candidate_gene_id"]) for r in entries}
        prediction_regions[key] = {r["region_id"] for r in entries}

    flow = Counter(total_source_records=len(source))
    scored = []
    excluded = []
    for number, row in enumerate(source, 1):
        pair_id = f"EPCB_{number:05d}"
        chromosome = row["chrom"]
        start, end = int(row["chromStart"]), int(row["chromEnd"])
        gene_external = row["measuredGeneEnsemblId"].split(".")[0]
        symbol = row["measuredGeneSymbol"]
        if gene_external in genes:
            mapping = "DIRECT_ENSEMBL_ID_GENCODE50"
            mapped = gene_external
        elif len(gene_symbols.get(symbol, set())) == 1:
            mapping = "UNIQUE_SYMBOL_TO_GENCODE50_FALLBACK"
            mapped = next(iter(gene_symbols[symbol]))
        elif len(gene_symbols.get(symbol, set())) > 1:
            mapping = "AMBIGUOUS_SYMBOL"
            mapped = "NA"
        else:
            mapping = "UNRESOLVED_IDENTIFIER"
            mapped = "NA"
        if mapping.startswith("DIRECT") and symbol not in genes[gene_external]:
            mapping = "DIRECT_ENSEMBL_ID_SYMBOL_MISMATCH"
        base = {
            "benchmark_pair_id": pair_id, "chromosome": chromosome, "start": start, "end": end,
            "external_gene_identifier": row["measuredGeneEnsemblId"], "external_gene_symbol": symbol,
            "mapped_gene_id_base": mapped, "mapping_status": mapping,
            "mapping_method": "EXTERNAL_ENSEMBL_BASE_ID_TO_GENCODE50",
            "external_label": "POSITIVE" if row["Regulated"] == "TRUE" else "NEGATIVE",
            "significant": row["Significant"], "regulated": row["Regulated"],
            "effect_size": row["EffectSize"], "adjusted_pvalue": row["pValueAdjusted"],
            "source_study": row["Dataset"], "source_reference": row["Reference"],
            "cell_type": row["CellType"], "genome_build": "GRCh38",
            "external_distance_to_tss": row["distanceToTSS"], "coordinate_system": "0-based-half-open",
        }
        reason = None
        if row["CellType"] == "K562": flow["k562_records"] += 1
        else: reason = "NON_K562"
        if reason is None and chromosome in ALLOWED and chromosome == row["chrTSS"]:
            flow["build_and_chromosome_compatible_records"] += 1
        elif reason is None: reason = "BUILD_OR_CHROMOSOME_UNRESOLVED"
        if reason is None and mapped != "NA": flow["mapped_gene_records"] += 1
        elif reason is None: reason = "GENE_ID_UNRESOLVED"
        if reason is None and row["ValidConnection"] == "TRUE": flow["valid_connection_records"] += 1
        elif reason is None: reason = "INVALID_OR_AMBIGUOUS_CONNECTION"
        promoter_overlap = False
        if reason is None:
            promoter_overlap = any(chromosome == pchrom and pend > start and pstart < end for pstart, pend, pchrom in promoters_by_gene[mapped])
            if promoter_overlap: reason = "PROMOTER_PROXIMAL_EXTERNAL_TEST"
            else: flow["distal_eligible_records"] += 1
        if reason:
            excluded.append({**base, "exclusion_reason": reason})
            continue

        region_hits = overlaps(region_index, region_max, chromosome, start, end)
        hit_ids = sorted({r["region_id"] for r in region_hits})
        midpoint = (start + end) // 2
        midpoint_ids = sorted({r["region_id"] for r in region_hits if int(r["start"]) <= midpoint < int(r["end"])})
        best_distance, nearest_genes = nearest_tss(tss_by_chr, chromosome, start, end)
        record = {**base, "external_positive": row["Regulated"] == "TRUE",
                  "overlapping_region_ids": ";".join(hit_ids) or "NA", "overlap_region_count": len(hit_ids),
                  "midpoint_region_ids": ";".join(midpoint_ids) or "NA", "nearest_tss_distance": best_distance,
                  "nearest_tss_gene_match": mapped in nearest_genes,
                  "colocated_region_overlap": any(region_id in colocated_atac for region_id in hit_ids)}
        for distance in distances:
            matching = sorted(rid for rid in hit_ids if (rid, mapped) in predictions[distance])
            midpoint_matching = sorted(rid for rid in midpoint_ids if (rid, mapped) in predictions[distance])
            record[f"candidate_{distance}"] = bool(matching)
            record[f"candidate_midpoint_{distance}"] = bool(midpoint_matching)
            record[f"candidate_region_ids_{distance}"] = ";".join(matching) or "NA"
            record[f"distance_only_{distance}"] = best_distance <= distance and mapped in nearest_genes
            record[f"activity_region_{distance}"] = bool(hit_ids)
            record[f"colocated_candidate_{distance}"] = any(rid in colocated_atac for rid in matching)
        scored.append(record)

    flow["final_scored_records"] = len(scored)
    flow["experimentally_labelled_positive_records"] = sum(r["external_positive"] for r in scored)
    flow["experimentally_labelled_negative_records"] = sum(not r["external_positive"] for r in scored)

    metric_rows = []
    for distance in distances:
        for method, subset in ((f"candidate_{distance}", "ALL_DISTAL_ELIGIBLE"),
                               (f"distance_only_{distance}", "ALL_DISTAL_ELIGIBLE"),
                               (f"activity_region_{distance}", "ALL_DISTAL_ELIGIBLE"),
                               (f"candidate_midpoint_{distance}", "MIDPOINT_MAPPING_SENSITIVITY"),
                               (f"colocated_candidate_{distance}", "ATAC_H3K27AC_COLOCATED_CANDIDATE")):
            metric_rows.append(metric_row(method, distance, subset, scored))
        for study in sorted({r["source_study"] for r in scored}):
            selected = [r for r in scored if r["source_study"] == study]
            metric_rows.append(metric_row(f"candidate_{distance}", distance, "BY_SOURCE_STUDY", selected, study))

    strata = []
    for low, high, label in DISTANCE_BINS:
        selected = [r for r in scored if low < abs(int(r["external_distance_to_tss"])) <= high or (low == 0 and abs(int(r["external_distance_to_tss"])) <= high)]
        if selected:
            strata.append(metric_row("candidate_100000", 100000, label, selected))

    failures = []
    for row in scored:
        predicted = row["candidate_100000"]
        if row["external_positive"] and predicted: classification, failure = "TRUE_POSITIVE", "RECOVERED_TESTED_POSITIVE"
        elif (not row["external_positive"]) and predicted: classification, failure = "FALSE_POSITIVE", "CANDIDATE_OVERLAPS_EXPERIMENTALLY_NEGATIVE_PAIR"
        elif (not row["external_positive"]) and (not predicted): classification, failure = "TRUE_NEGATIVE", "NOT_PREDICTED_TESTED_NEGATIVE"
        elif not row["overlap_region_count"]: classification, failure = "FALSE_NEGATIVE", "ELEMENT_NOT_REPRESENTED_BY_ACTIVITY_REGION"
        elif int(row["external_distance_to_tss"]) > 100000: classification, failure = "FALSE_NEGATIVE", "TESTED_GENE_BEYOND_100KB"
        elif not row["nearest_tss_gene_match"]: classification, failure = "FALSE_NEGATIVE", "TESTED_GENE_NOT_NEAREST_TSS"
        else: classification, failure = "FALSE_NEGATIVE", "ACTIVITY_REGION_PRESENT_LINK_NOT_RECOVERED"
        failures.append({"benchmark_pair_id": row["benchmark_pair_id"], "classification": classification,
                         "failure_class": failure, "source_study": row["source_study"],
                         "external_label": row["external_label"], "predicted_positive": str(predicted).lower(),
                         "overlapping_region_ids": row["overlapping_region_ids"],
                         "mapped_gene_id_base": row["mapped_gene_id_base"]})

    public_fields = [k for k in scored[0] if k != "external_positive"] if scored else []
    for row in scored:
        for key, value in list(row.items()):
            if isinstance(value, bool): row[key] = str(value).lower()
    write_tsv(out / "crispr_benchmark_pairs.tsv", public_fields, scored)
    write_tsv(out / "crispr_benchmark_exclusions.tsv", list(excluded[0]) if excluded else ["benchmark_pair_id", "exclusion_reason"], excluded)
    write_tsv(out / "crispr_benchmark_denominators.tsv", ["stage", "n_records"], [{"stage": k, "n_records": v} for k, v in sorted(flow.items())])
    metric_fields = list(metric_rows[0])
    write_tsv(out / "crispr_benchmark_metrics.tsv", metric_fields, metric_rows)
    write_tsv(out / "distance_stratified_metrics.tsv", metric_fields, strata)
    write_tsv(out / "crispr_failure_classes.tsv", list(failures[0]), failures)
    counts = Counter(row["failure_class"] for row in failures)
    write_tsv(out / "crispr_failure_summary.tsv", ["failure_class", "n_records"], [{"failure_class": k, "n_records": v} for k, v in sorted(counts.items())])
    print(f"source={len(source)} scored={len(scored)} positives={flow['experimentally_labelled_positive_records']} negatives={flow['experimentally_labelled_negative_records']} excluded={len(excluded)}")


if __name__ == "__main__":
    main()
