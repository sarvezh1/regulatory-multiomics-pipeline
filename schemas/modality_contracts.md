# Standardized modality contracts

All tabular outputs are UTF-8 TSV with a header, explicit `NA`, stable sort order, and a sidecar provenance record. Genomic interval tables declare build and coordinate system. Statistical tables include contrast and analysis state.

The production annotation is GENCODE 50, based on GRCh38.p14, while the production coordinate representation is a separately pinned GRCh38 primary-assembly FASTA. Default integration is limited to primary chromosomes `chr1`–`chr22`, `chrX`, and `chrY`. Patches, alternate loci, haplotypes, unlocalized/unplaced sequences, and mitochondrial sequence are excluded from default regulatory integration.

Every externally processed input must carry the following provenance fields, using `UNKNOWN` when absent from the source and `UNRESOLVED` when ambiguous or unverified:

- `source_genome_build`
- `source_annotation`
- `source_processing_pipeline`
- `source_file_accession`
- `source_file_checksum`
- `coordinate_normalization_status`
- `gene_id_mapping_status`

Externally processed files are never relabelled as GENCODE 50 merely because GENCODE 50 is the production annotation. Original identifiers and, where applicable, original coordinates remain recoverable after harmonization.

## RNA

Required logical outputs: sample metadata; raw integer gene counts; normalized descriptive abundance; differential-expression statistics; original, stable, versioned, and harmonized gene identifiers as applicable; source annotation and mapping status; sample and library QC. Normalized abundance is never substituted for raw counts in count-model differential analysis.

## ATAC

Required logical outputs: frozen consensus genomic regions with stable coordinate-derived identifiers and contributing-peak provenance; one integer paired-fragment count per biological sample and region; explicitly named descriptive normalization; complete differential-accessibility statistics when valid; sample peaks; library/sample metadata; peak, alignment, complexity, fragment-length, TSS, FRiP, and replicate QC; component status; and reference/tool provenance. Raw technical libraries remain distinct through QC and are combined before sample-level peaks and counts. Source and normalized genome build, coordinate system, and normalization status remain explicit.

## Methylation

The implemented WGBS CpG contract uses deterministic `cpg_id` values derived from `chromosome:start-end`, 0-based half-open two-base CpG intervals, methylated and unmethylated integer counts, coverage, methylation fraction, source-strand provenance, measurement state, and original coordinate provenance. Coverage always equals methylated plus unmethylated counts. Complementary strands are merged by summing counts before computing the fraction. Calls that do not resolve to a reference CpG are retained in an explicit exclusion table rather than entering the canonical contract.

Technical-library counts are summed to biological samples. The cross-sample long table preserves uncovered loci as `MISSING` with `NA` counts and fractions; zero coverage is never represented as zero methylation. Count-preserving processed input converges on this contract. Percentage-only input is descriptive and cannot enter dmrseq.

Region-level differential output records deterministic DMR identifiers, 0-based half-open coordinates, CpG count, numerator and denominator means, fraction-scale and percentage-point deltas, statistics, p/q values, contrast, and method. Positive delta means greater methylation in the named numerator.

## H3K27ac

Required logical outputs: enriched regions; per-sample integer counts and/or clearly labeled descriptive signal; differential-enrichment statistics when supported; input-control identifiers and provenance; source and normalized genomic coordinates; FRiP/cross-correlation/library QC where available. Signal-only tracks are valid descriptive inputs but do not satisfy a count-based differential contract. Missing input-control identity blocks raw H3K27ac validation but does not block processed descriptive integration or unrelated modality development.

## Integrated

The implemented contract is relational. `promoter_evidence.tsv` is the central promoter summary, while ATAC–promoter, H3K27ac–promoter, CpG–promoter, and DMR–promoter tables preserve every detailed one-to-many relationship. Companion tables preserve gene-ID mappings, promoter methylation by sample, gene summaries, ATAC/H3K27ac co-localizations, nearest-TSS candidates, and external enhancer–gene evidence.

Every promoter summary records explicit matching class, assay availability, assay test state, effects and assay-specific statistics where valid, literal regulatory-pattern state, compatible contrast, component status, and provenance identifier. Original and base Ensembl identifiers coexist. Missing modalities, untested features, nonsignificant observations, and observed zeros remain distinct. Empty valid results and unrun analyses are represented by state, never fabricated rows or zeros. Exact fields are recorded in `schemas/integration_contract.tsv`; scientific semantics are defined in `docs/integration_branch.md` and the frozen `docs/integration_specification.md`.
