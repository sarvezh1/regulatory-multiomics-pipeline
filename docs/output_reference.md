# Output reference

All contracts use tab-separated text unless noted. Companion provenance records the schema version, parameters, references, tools and component statuses. Missing values remain `NA`, `UNKNOWN`, `UNRESOLVED`, `NOT_AVAILABLE` or `NOT_TESTED` according to the contract; they are never silently replaced by biological zero.

## External-validation outputs

External-validation outputs are companion evidence rather than production assay contracts. `ccre_concordance_summary.tsv` and `ccre_region_overlaps.tsv` retain source cCRE categories and strict positive-base overlaps. `crispr_benchmark_pairs.tsv` has one distal-eligible experimentally tested element–gene pair per row, keyed by `benchmark_pair_id`, with original and mapped identifiers, source study, external label, region overlaps, and predictions at frozen distances. `crispr_benchmark_metrics.tsv` records method, distance, subset, denominator, class counts, confusion cells and binary metrics. `crispr_failure_classes.tsv` preserves row-level TP/FP/TN/FN classification and failure class. `biological_plausibility_summary.tsv` and `biological_plausibility_correlations.tsv` are descriptive `REFERENCE_COMPARABLE` summaries.

All external-validation genomic intervals are GRCh38 0-based half-open. `external_validation_status.tsv` retains validation-component status separately from upstream production statuses. `external_validation_provenance.json` records the schema version, frozen rules, independence classes, checksums and optional-analysis states. Full definitions and claim boundaries are in [external validation](external_validation.md).

## RNA

`rna/counts/rna.raw_counts.tsv` is gene-by-biological-sample integer count data keyed by preserved original and base Ensembl gene identifiers. Descriptive normalized values are not inferential inputs. Differential output retains the complete tested-gene table, positive log2 fold change meaning greater expression in the named numerator, and preserves missing statistics. See `schemas/modality_contracts.md` and `docs/rna_branch.md`.

## ATAC

`atac/counts/atac.raw_counts.tsv` contains one row per deterministic 0-based half-open consensus region and one integer fragment count per biological sample. Detailed peaks, sample metadata, QC and optional complete differential tables are companion outputs. Positive differential effect means greater accessibility in the numerator. See `docs/atac_branch.md`.

## WGBS

`wgbs/counts/wgbs.cpg.long.tsv` contains canonical CpG loci, methylated and unmethylated counts, total coverage and fractions. Loci use 0-based half-open intervals; coverage equals methylated plus unmethylated counts. No coverage remains missing, not zero methylation. Region inference uses dmrseq; positive delta means greater methylation in the numerator. See `docs/wgbs_branch.md`.

## H3K27ac

`h3k27ac/counts/h3k27ac.raw_counts.tsv` contains deterministic consensus regions and integer treatment evidence per biological sample. Control relationships, sample peaks, signal, QC and optional complete differential tables remain traceable. Positive differential effect means greater H3K27ac enrichment in the numerator. See `docs/h3k27ac_branch.md`.

## Integration

`integration/evidence/promoter_evidence.tsv` is keyed by promoter identity and preserves assay availability, tested state, effects and literal evidence patterns. Detailed ATAC-promoter, H3K27ac-promoter, CpG-promoter and DMR-promoter tables preserve one-to-many relationships. Gene summaries retain promoter heterogeneity. Distal nearest-TSS links are predictions, not validated targets. See `schemas/integration_contracts.md` and `docs/integration_branch.md`.

## Status and provenance

Analysis status files are component scoped. `COMPLETE_EMPTY_RESULT` is a successful analysis with no result; `NOT_RUN_*` and `NOT_AVAILABLE_INPUT` are not tool crashes; `FAILED` denotes an actual failure. Unified manifest formats and primary keys are listed in `schemas/run_provenance_contract.tsv`.
