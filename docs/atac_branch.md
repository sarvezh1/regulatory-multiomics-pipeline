# ATAC-seq branch

Author: Sarvesh Galgale

## Scope and inputs

The v0.1 raw route supports paired-end human bulk ATAC-seq only. Each row is a technical library (`library_id`); `sample_id` is the biological sample. Technical libraries receive independent raw QC and are deterministically combined before sample-level peak calling, consensus counting, and differential inference. They never increase biological n. Accepted assay labels are exactly `ATAC`, `ATAC-seq`, `atac`, and `atac-seq`, normalized internally to `ATAC`.

The processed route accepts one tab-separated common-region integer count matrix with `region_id`, `chromosome`, `start`, `end`, and one column per biological sample. Rows must provide explicit source build, pipeline, accession, checksum, and coordinate-normalization provenance or an explicit missing-provenance state. Peak-only input is deferred: peak scores are not pseudo-counts and cannot support quantitative differential analysis.

## Reference and coordinates

Production requires a separately pinned GRCh38 primary-assembly FASTA, matching GENCODE 50 GTF, and GRCh38 ENCODE blacklist `ENCFF356LFX`. Default nuclear analysis uses chromosomes 1–22, X, and Y. `chrM` is counted for QC and excluded before nuclear peak analysis. Patches, alternate loci, haplotypes, and unplaced/unlocalized sequences are outside the default integration space.

The TSS BED is derived from transcript records in the supplied GTF. Identical chromosome/strand/TSS positions are deduplicated, stored as one-base zero-based half-open intervals, and retain transcript/gene provenance. TSS QC windows are independent of the regulatory promoter definition.

The compact chromosome 22 reference is execution-only. Its Bowtie2 index, TSS BED, and blacklist subset are derived deterministically during validation; it is not a production or biological benchmark.

## Raw route

Libraries pass through raw FastQC, conservative fastp preprocessing, post-trim FastQC, and Bowtie2 with `--very-sensitive -X 2000 --no-mixed --no-discordant`. The filter sequence is: map; name-sort; `samtools fixmate -m`; coordinate-sort; duplicate marking; primary/proper-pair and MAPQ >=30 selection; mitochondrial exclusion; duplicate exclusion; then fragment-level blacklist exclusion. Counts for unmapped, mapped, proper, high-quality, mitochondrial, duplicate, blacklist, and usable fragments remain machine-readable.

Canonical BAM coordinates are never globally shifted. A +4/-5 insertion representation is derived only for TSS QC. Duplicate marking is library-local. Sample BAMs and fragments combine technical libraries only within `sample_id`.

MACS3 uses `-f BAMPE -q 0.01 --keep-dup all --call-summits` on the already deduplicated, mitochondrial-free, blacklist-free sample BAM. It uses observed fragments, narrow peaks, no single-end shift, no conventional input control, and an explicit effective genome size. CPM bigWigs represent extended paired-fragment coverage and are descriptive, never the DA input.

## Consensus and counts

Valid sample narrow peaks are sorted deterministically. Overlap-connected intervals are merged; book-ended intervals remain separate. Each consensus ID is `ATAC_` plus the first 16 hexadecimal characters of SHA-256 over `chromosome:start-end`. A provenance table maps every region to its contributing sample peaks. One paired fragment is counted once for each consensus region it overlaps; technical libraries contribute to their biological sample column.

## QC and differential analysis

Outputs include FastQC/fastp/Bowtie2/samtools metrics, mapped/proper/mitochondrial/duplicate fractions, usable fragments, fragment-length distribution, nucleosome-free and mono-/multinucleosomal counts, TSS enrichment, FRiP, peak count, NRF/PBC1/PBC2, sample-peak Jaccard, count correlations, and PCA where estimable. ataqv JSON remains available alongside the deterministic tables. MultiQC parses the supported logs locally; publication retains a deterministic static summary and curated parsed TSV tables. No single QC threshold automatically excludes a sample.

DA is opt-in with `--run_da true`, a safe formula such as `~ condition` or an estimable `~ batch + condition`, and `--contrast variable,numerator,denominator`. Arbitrary R expressions are rejected. Each compared group needs at least two independent biological samples; n=2 is permitted but weak and n>=3 is preferred. Positive log2 fold change means greater accessibility in the named numerator. Complete tested-region tables preserve `NA`.

Statuses are component scoped: `atac_processing`, `atac_peak_calling`, `atac_count_contract`, and `atac_differential`. Zero valid peaks are `COMPLETE_EMPTY_RESULT`; disabled or scientifically invalid inference uses the appropriate `NOT_RUN_*` state; an actual tool failure remains a workflow failure.

## Limitations

The branch does not reproduce the full ENCODE uniform ATAC pipeline. Raw single-end ATAC, peak-only quantitative input, IDR/pseudoreplicates, motifs, footprinting, ATAC annotation, and cross-modality integration are not implemented. Compact-reference QC values are execution diagnostics only.
