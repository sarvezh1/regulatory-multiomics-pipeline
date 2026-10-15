# H3K27ac ChIP-seq branch

Author: Sarvesh Galgale

## Scope

The v0.1 histone branch supports human bulk H3K27ac ChIP-seq only. H3K27ac is treated as a narrow mark. Raw single-end and paired-end routes are implemented and validated. CUT&Tag, other histone targets, transcription-factor ChIP, control-free peak calling, IDR, motif analysis, and cross-modality integration are outside this block.

Raw treatment rows use assay `H3K27ac_ChIP` and target `H3K27ac`. Input rows use assay `ChIP_input` and target `CONTROL`. Every raw treatment record must name a compatible input record through `control_sample_id`; controls are never inferred from filenames, row order, or accession similarity. Compatibility validation requires matching condition, biosample label, and layout. `sample_id` denotes biological material and `library_id` denotes a technical library.

## Raw route

FastQC 0.12.1 assesses raw and trimmed reads. fastp 0.23.4 applies its non-aggressive default quality and adapter handling to treatment and control alike. Bowtie2 2.5.5 uses `--very-sensitive` for single-end libraries and `--very-sensitive --no-mixed --no-discordant -X 1000` for paired-end libraries. The paired fragment maximum is ChIP-specific and does not reuse the ATAC 2,000-bp policy.

Each technical library is duplicate-marked independently with samtools 1.20. Paired-end preparation is name sort, `fixmate -m`, coordinate sort, and `markdup`; single-end preparation is coordinate sort and `markdup`. Primary mapped records with MAPQ at least 30 are retained, paired-end evidence additionally requires proper pairing, duplicates are excluded, chrM is counted then excluded, and GRCh38 blacklist overlaps are removed. Input and treatment use the same sequence. NRF, PBC1, and PBC2 are descriptive and do not impose automatic pass/fail decisions.

Technical libraries are combined only after library-level filtering. Peak calling occurs at the biological-sample level with its explicit matched input.

## Peak calling and signal

MACS3 3.0.5 uses a narrow-peak q-value threshold of 0.01, `--keep-dup all` after upstream duplicate removal, `--call-summits`, `-B`, and `--SPMR`. Single-end data use `-f BAM` and normal model construction by default. Paired-end data use `-f BAMPE` and observed fragments; no single-end shift or extension is applied. Effective genome size is explicit and must match the reference scope.

The optional `chip_se_validation_extsize` parameter defaults to zero and therefore leaves the production model enabled. A nonzero value is reserved for a labelled compact execution fixture when an independently documented fragment extension is available; it activates `--nomodel --extsize`. It must not be tuned to create desired peaks.

Outputs include narrowPeak, summits, MACS3 pileup and control-lambda bedGraphs, MACS3 `bdgcmp -m FE` fold enrichment, treatment and input CPM bigWigs, and a deepTools CPM-normalized treatment/input log2-ratio bigWig. Signal tracks are descriptive and never replace integer counts for inference.

FRiP is sample peak-overlapping evidence divided by usable nonduplicate nuclear reads for single-end data or fragments for paired-end data after MAPQ, mitochondrial, duplicate, and blacklist filtering. NSC and RSC are explicitly `NOT_IMPLEMENTED`; a reproducible cross-correlation dependency was not added merely to approximate those measures.

## Consensus regions and counts

Sample narrow peaks are sorted and merged through connected strict overlaps. Book-ended intervals remain separate. Every region retains contributing sample-peak provenance and receives a coordinate-derived `H3K27AC_` identifier. Single-end reads are counted once; paired-end fragments are counted once. Technical libraries contribute to one biological-sample column. Normalized bigWig values are not count inputs.

## Processed route and differential analysis

The quantitative processed route accepts only `h3k27ac_region_count_matrix_tsv`: valid 0-based half-open coordinates, stable region identifiers, non-negative integer treatment/control counts, biological-sample metadata, target/control identity, GRCh38 compatibility, and source provenance. It does not infer counts from peak scores, q-values, fold enrichment, or bigWigs. Peak-only quantitative inference is therefore refused; descriptive peak-only adaptation is deferred.

Differential H3K27ac enrichment is opt-in with `--run_chip_diff true`, a safe formula such as `~ condition` or an estimable `~ batch + condition`, and `--contrast variable,numerator,denominator`. DiffBind 3.20.0 provides the sample/region framework and DESeq2 1.50.2 the validated integer-count model. Matched controls define raw peaks and remain explicit provenance; counts are not manually control-subtracted. At least two independent biological samples are required per group, with three preferred. Positive log2 fold change means greater H3K27ac enrichment in the named numerator. Complete tested-region tables preserve NA values.

GSE152146 H3K27ac remains descriptive because relevant conditions lack independent replication. This implementation does not change that frozen decision.

## Statuses and limitations

Statuses are component-scoped: `h3k27ac_processing`, `h3k27ac_peak_calling`, `h3k27ac_count_contract`, and `h3k27ac_differential`. Missing control is `NOT_AVAILABLE_INPUT`; a legitimate empty peak set is `COMPLETE_EMPTY_RESULT`; a tool crash is `FAILED`. Disabling differential analysis does not change a successful processing status.

The branch does not reproduce the complete ENCODE uniform ChIP-seq/IDR pipeline. Compact chr22 diagnostics are not genome-wide QC or biological benchmarks. Replicate overlap/correlation is supported structurally by common-region counts but no automatic sample exclusion threshold is imposed.
