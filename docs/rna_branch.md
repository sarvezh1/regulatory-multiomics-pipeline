# RNA-seq branch

Author: Sarvesh Galgale

## Scope and architecture

The RNA branch implements only bulk short-read RNA-seq. `main.nf` invokes the DSL2 RNA workflow in `workflows/rna`; individual processes live in `modules/local`. An RNA samplesheet may contain raw libraries, processed gene-count libraries, or both. Any assay other than the exact accepted values `RNA`, `RNA-seq`, `rna`, and `rna-seq` is rejected by the RNA branch rather than ignored.

The production reference contract is Homo sapiens, a separately pinned GRCh38 primary-assembly FASTA, and GENCODE release 50 annotation (GRCh38.p14). Production reference validation requires exactly `chr1`-`chr22`, `chrX`, `chrY`, and `chrM`; it checks that every GTF sequence name exists in the FASTA and rejects alternate, patch, haplotype, unplaced, and unlocalized names. The annotation release and FASTA representation remain distinct provenance fields. The chr22 reference under `tests/integration/reference` is execution-only and is not a production reference or biological benchmark.

## Samples and libraries

`sample_id` is the biological sample identity. The RNA branch requires `library_id` because one biological sample can have multiple technical libraries without acquiring additional biological replication. Each library runs through its own raw QC and processing path. Integer gene counts from libraries sharing a `sample_id` are summed deterministically before DESeq2. The technical-replicate manifest records every library and the combination rule.

`strandedness` is required for every RNA row and accepts only `unstranded`, `forward`, or `reverse`. It maps to featureCounts `-s 0`, `-s 1`, or `-s 2`. Strandedness is never guessed. Raw single-end rows require `read1` and forbid `read2`; raw paired-end rows require both. Processed rows use `library_layout=not_applicable`, provide no reads, and use `processed_input_type=gene_count_tsv`.

## Raw route

The raw route executes:

1. FastQC on source reads;
2. fastp with non-aggressive defaults and JSON/HTML metrics;
3. FastQC on retained reads;
4. deterministic STAR index generation from the supplied FASTA/GTF;
5. STAR splice-aware alignment to coordinate-sorted BAM;
6. samtools BAM indexing, flagstat, and stats;
7. featureCounts over `exon`, grouped by `gene_id`, with explicit layout and strandedness;
8. adaptation to the RNA count contract;
9. MultiQC parsing, deterministic static QC summary, and provenance.

STAR standard multimapping behavior is retained. On Windows-hosted Docker/WSL execution, STAR's FIFO-only temporary directory is placed in container-local `/tmp`; scientific outputs remain in the Nextflow work directory. The workflow publishes retained deliverables by copy because Windows bind mounts do not provide a uniformly reliable link mode. Temporary work products remain resumable and may be removed after successful review.

## Processed route and identifiers

The canonical external format is a UTF-8 tab-separated file containing `gene_id_original` and non-negative integer `raw_count`, with optional `gene_name` and `gene_id_base`. The adapter derives `gene_id_base` only by removing a terminal version component when the source does not provide it. It performs no gene-symbol mapping. The original identifier remains in every downstream table.

Processed rows must explicitly provide source build, annotation, processing pipeline, file accession, checksum, coordinate-normalization status, and gene-ID mapping status. `UNKNOWN` means absent from the source; `UNRESOLVED` means present but not confidently resolved. External counts are never relabelled as GENCODE 50.

## Differential expression

DESeq2 uses the combined integer count matrix. Descriptive output is named `rna.normalized_counts.tsv` and contains DESeq2 size-factor-normalized counts; it is not TPM. Inference requires `--run_de true` and `--contrast variable,test_level,reference_level`. `--design_formula` accepts only additive column names, for example `~ condition` or `~ batch + condition`; arbitrary R expressions are rejected.

Each compared level requires at least two independent biological samples and the model matrix must be full rank. Two per group permits fitting but has limited power; three or more is preferred. Complete tested-gene results retain `NA` and include `gene_id_original`, `gene_id_base`, `gene_name`, `baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, and `padj`. A valid analysis with no genes below FDR 0.05 receives `COMPLETE_EMPTY_RESULT`, not an unrun state.

## Outputs and limitations

Outputs are organized under `rna/qc`, `rna/trimmed`, `rna/alignment`, `rna/counts`, `rna/differential`, `rna/multiqc`, and `rna/provenance`. MultiQC parses supported logs locally; retained publication consists of a static HTML summary and curated parsed TSV tables. The provenance manifest removes machine-specific input paths, records basenames and SHA-256 checksums, reference checksums/source URLs, samples/design/contrast, exact executed versions, pinned containers, run time, and analysis states.

The RNA branch does not implement transcript-level quantification, isoform inference, ATAC-seq, WGBS, H3K27ac, regulatory integration, or publication-scale reporting. A full production GRCh38 STAR index was not built locally.
