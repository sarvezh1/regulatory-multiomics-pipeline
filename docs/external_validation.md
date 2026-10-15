# External validation

External validation evaluates the frozen v0.1 regulatory-linking behavior; it does not tune or replace it. Production candidates are generated before external labels are loaded, and their checksums are recorded in `metadata/external_validation_prediction_freeze.tsv`. The primary 100 kb nearest-TSS rule and the predefined 50 kb and 250 kb sensitivities remain unchanged regardless of benchmark performance.

## Sources and independence

Exact source versions, accessions, retrieval dates, checksums, usage status and validation roles are recorded in `metadata/external_validation_sources.tsv`.

- ENCODE K562 ATAC `ENCFF057UYP` and H3K27ac `ENCFF038DDS` released GRCh38 replicated-peak files form a validation-only full-genome regulatory-region set. They are independent experiments and remain `REFERENCE_COMPARABLE`, not matched aliquots.
- ENCODE Encyclopedia v4 K562 cCRE file `ENCFF286VQG` from `ENCSR913HQX` is used for annotation concordance. Its classification is `ANNOTATION_CONCORDANCE_NONINDEPENDENT` because the catalogue incorporates chromatin evidence that can overlap the source assays.
- `EPCrisprBenchmark_combined_data.training_K562.GRCh38.tsv.gz` from EngreitzLab/CRISPR_comparison commit `c075c1317a5d8bb47ad2fab33c8eb9d4bce99aac` is the independent functional benchmark. Source-study labels are retained.
- GENCODE release 50 annotation and its separately identified GRCh38.p14 primary-assembly FASTA provide the full-genome TSS and chromosome-boundary references. Annotation release and sequence representation are not treated as interchangeable.

Downloaded source resources are local validation inputs and are excluded from future repository inclusion. Fetch locations and checksums, rather than third-party data copies, are public project records.

## Frozen prediction construction

The validated integration rules are replayed on full GRCh38 primary chromosomes 1–22, X and Y. All distinct GENCODE v50 transcript TSS positions are retained, with exact same-gene/same-strand duplicates collapsed. Promoters use the strand-aware -2000/+500 definition and strict positive-base overlap. Released ATAC regions are the distal anchors; overlapping H3K27ac is retained as co-evidence. The nearest TSS is linked only within the configured distance, and exact ties are retained.

The streaming full-reference builder was verified field-for-field against the integration chr22 reference output for all 8,676 TSS records and all three sensitivity promoter resources. It was used because loading the full annotation through the containerized reference builder exceeded the practical local validation runtime and memory envelope; this changed implementation scale, not scientific rules.

## cCRE annotation concordance

The current cCRE BED schema was audited before analysis. Its exact category vocabulary is retained: `PLS`, `pELS`, `dELS`, `CA-CTCF`, `CA-H3K4me3`, `CA-TF`, `CA-only`, and `Low-DNase`. No categories are renamed as active enhancers. Strict interval intersections are reported for promoter ATAC, promoter H3K27ac, distal ATAC-only, distal H3K27ac-only, and distal ATAC/H3K27ac-colocated regions. No significance test or combined validation score is assigned.

## CRISPR functional benchmark

The source defines `Regulated=TRUE` as an adjusted p-value below 0.05 with a negative target-gene expression effect. `Regulated=FALSE` is the negative class among source-defined valid connections. Ambiguous or invalid labels would be excluded; the pinned file marks every row as a valid connection.

External base Ensembl identifiers are joined directly to GENCODE v50. A unique current GENCODE symbol is used only as a recorded fallback when the external Ensembl identifier is absent from release 50; ambiguous and unresolved mappings are excluded. Original identifiers and symbols remain in the pair table.

The primary scoring subset is K562, GRCh38, primary chromosomes, mapped genes, valid source connections, and distal tests. An external element with any positive-base overlap to a frozen promoter of its tested gene is classified `PROMOTER_PROXIMAL_EXTERNAL_TEST` and excluded from distal scoring. No liftover is performed.

A pair is predicted positive only when the tested element overlaps a frozen K562 distal region and that same region has a nearest-TSS candidate link to the tested gene within the tested distance. The primary mapping uses any positive-base intersection; the predefined sensitivity requires the external element midpoint to fall in the predicted region.

Reported binary measures are TP, FP, TN, FN, precision, recall, specificity, F1, balanced accuracy and Matthews correlation coefficient. Fisher's exact odds ratio and two-sided p-value are supplementary; the approximate 95% odds-ratio interval is reported only when all cells are nonzero. No ROC-AUC or PR-AUC is assigned to the binary predictor.

Two validation-only baselines separate region detection and gene assignment: distance-only nearest TSS ignores chromatin evidence, while activity-region overlap ignores the tested gene. Neither baseline changes the production method. Failure tables distinguish missing activity representation, tested genes that are not nearest, tested genes beyond 100 kb, and residual link non-recovery.

## Biological plausibility

The existing K562 integration remains `REFERENCE_COMPARABLE`. Genes are split using frozen inclusive 25th and 75th percentile expression boundaries into bottom, middle and top strata. Ties at zero expression make stratum sizes unequal and are retained deterministically. Promoter ATAC, promoter H3K27ac and coverage-weighted promoter methylation summaries are descriptive. Spearman rho is reported across genes with complete observations, without a correlation p-value and without matched-sample or causal claims.

## Deferred optional comparisons

No ABC model is run and no Hi-C data are downloaded. The optional ABC comparator is `NOT_RUN_DISABLED`. The optional GSE152146 RNA+ATAC differential audit is `NOT_AVAILABLE_INPUT` with reason `DEFERRED_NO_VALID_COMPACT_QUANTITATIVE_INPUT`; normalized signal products are not substituted for count-preserving inputs.

Machine-readable outputs and exact headline sources are documented in the [external validation record](validation/external.md). The CRISPR benchmark covers experimentally tested loci and genes, not the whole genome.
