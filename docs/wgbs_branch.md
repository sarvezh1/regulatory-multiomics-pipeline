# Conventional WGBS branch

## Implemented scope

The v0.1 raw route supports conventional Illumina paired-end whole-genome bisulfite sequencing only. It rejects RRBS, enzymatic methyl-seq, oxidative bisulfite sequencing, PBAT, NOMe-seq, targeted bisulfite assays, and unspecified protocols. `bisulfite_directionality` is mandatory for raw data and accepts only `directional` or `non_directional`; it is never inferred from filenames or read orientation.

`sample_id` denotes a biological sample and `library_id` a technical library. Libraries receive independent QC, trimming, alignment, deduplication, and extraction. CpG counts are aggregated only after per-library deduplication. Technical libraries never increase biological replication.

## Raw route

The executed path is:

1. FastQC 0.12.1 on both raw mates.
2. Trim Galore 2.3.0 with paired mode, quality 20, minimum retained length 20, and no RRBS-specific or arbitrary fixed-base trimming.
3. Bismark 3.1.0 reference preparation with Bowtie2, using the supplied FASTA.
4. Bismark faithful/default paired short-read alignment. `--non_directional` is added only for rows explicitly declared non-directional.
5. Bismark paired deduplication per technical library before any aggregation.
6. Paired methylation extraction from deduplicated alignments with comprehensive context output, count output, M-bias output, and overlapping-mate suppression. The four configurable M-bias ignore values default to zero.
7. Reference-aware canonical CpG conversion, technical-library aggregation, cross-sample universe construction, QC, reporting, provenance, and optional dmrseq.

The default route does not produce a full all-cytosine report. CpG is the production contract; CHG and CHH are retained as summary QC. A low alignment rate is a QC observation, not a process failure.

## CpG and coordinate contract

Bismark coverage positions are interpreted as 1-based closed cytosine coordinates. Each call must resolve against the supplied FASTA to either the top-strand C or complementary bottom-strand G of a reference CpG. The canonical locus is a 0-based half-open two-base interval. Complementary evidence is merged by summing methylated and unmethylated counts before computing methylation fraction.

Polymorphic or otherwise non-reference CpG calls can occur in genuine reads. They are excluded from the reference CpG contract and retained with counts and rejection reason in `<library_id>.nonreference_cpg.tsv`. Impossible processed-input coordinates are rejected.

For every observed CpG:

`coverage = methylated_count + unmethylated_count`

Methylation fraction is `NA` when coverage is zero. In cross-sample tables, an uncovered CpG is `MISSING` with `NA` counts and fraction; it is not zero methylation. The default inferential coverage threshold is five reads, required in at least two independent samples in each compared group.

## Processed route

`cpg_count_tsv` requires chromosome, coordinates, coordinate convention, methylated count, and either unmethylated count or total coverage. If a fraction is supplied it must agree with the counts. Counts must be nonnegative, coordinates valid, and source build/provenance explicit.

`cpg_percentage_tsv` is accepted for descriptive output only. Counts are never reconstructed from percentages, and differential methylation reports `NOT_AVAILABLE_INPUT`. Count and percentage input tiers cannot be mixed in one run.

## Differential methylation

dmrseq 1.30.0 is the only production inferential method. It runs only when `--run_dm true`, an explicit safe formula is supplied, and `--contrast variable,numerator,denominator` names valid levels. Supported formulas include `~ condition` and estimable additive covariate designs such as `~ batch + condition`. Strings are parsed as variable names and are not evaluated as executable code.

Each group requires at least two independent biological samples; three or more are preferable. The design must be full rank and count-preserving CpGs must pass the frozen coverage rule. Missing measurements enter the internal `bsseq` object as zero coverage while remaining `MISSING` in the public contract. Positive `delta_methylation` means greater methylation in the numerator; values are fractions in [-1,1], with a separate percentage-point column.

dmrseq discovers candidate regions and does not use an ATAC-like fixed region universe. The complete returned table is retained. A valid table with no q-value below a reporting threshold remains `COMPLETE`; no candidate regions or no eligible CpGs use `COMPLETE_EMPTY_RESULT`.

## QC and conversion semantics

Machine-readable QC includes input/retained pairs, unique alignments, mapping efficiency, duplicates, CpG loci and coverage, excluded non-reference CpG calls, global CpG methylation, CHG/CHH methylation, M-bias status, and conversion-control status. A real named unmethylated control contig is required for conversion-efficiency reporting. Without one, observed non-CpG methylation is reported separately and conversion efficiency is `NOT_AVAILABLE_INPUT`.

## Statuses and limitations

Statuses are component scoped: `wgbs_processing`, `wgbs_methylation_extraction`, `wgbs_count_contract`, and `wgbs_differential`. Zero aligned reads can yield successful processing with `COMPLETE_EMPTY_RESULT` extraction. Tool crashes remain failures.

Raw single-end WGBS, non-conventional methylation protocols, CpG-level inferential DMLs, full production-reference execution, H3K27ac, and cross-modality integration are not implemented.

## Primary references

- [Bismark documentation](https://felixkrueger.github.io/Bismark/)
- [Trim Galore releases](https://github.com/FelixKrueger/TrimGalore/releases)
- [dmrseq 1.30.0](https://bioconductor.org/packages/3.22/bioc/html/dmrseq.html)
- [bsseq 1.46.0](https://bioconductor.org/packages/3.22/bioc/html/bsseq.html)
