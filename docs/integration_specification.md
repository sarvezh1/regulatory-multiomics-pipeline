# Integration specification

Status: frozen for v0.1

## Coordinate and identifier policy

- GENCODE release 50 is the production annotation and is based on GRCh38.p14. The annotation release is distinct from the production FASTA representation.
- The production FASTA is an independently pinned GRCh38 primary-assembly representation with checksum, sequence dictionary, and chromosome sizes. Default integration is restricted to `chr1`–`chr22`, `chrX`, and `chrY`. Assembly patches, alternate loci, haplotypes, unlocalized sequences, and unplaced sequences are excluded. `chrM` may contribute QC metrics but is excluded from default nuclear regulatory integration.
- A dataset label of `hg38` or `GRCh38` is insufficient to establish exact contig compatibility or GENCODE 50 annotation. Compatibility is verified from source provenance and sequence/contig metadata.
- Internal intervals use 0-based, half-open coordinates: `[start, end)`. Imported 1-based formats are converted once at a named boundary and retain the original coordinates, source build, conversion method, and conversion status in provenance.
- Sequence aliases are normalized only through a versioned mapping. Coordinate normalization never implies annotation or identifier normalization.
- Stable Ensembl gene IDs without version suffixes are primary keys. The original versioned ID is retained. Gene symbols are annotations, never join keys.
- Externally processed expression identifiers are never relabelled as GENCODE 50 identifiers. Any mapping to GENCODE 50 creates a new harmonized identifier field, records the mapping source/status, and preserves the original identifier and source annotation so the transformation is reversible.
- Each distinct annotated TSS is retained initially. The canonical promoter-level view uses one deterministic representative TSS per gene: GENCODE canonical transcript when present, then MANE Select, then longest protein-coding transcript, then lowest genomic TSS coordinate as a final tie-breaker. A sensitivity view uses all unique TSSs.

## Promoters and overlaps

- Default promoter: strand-aware `TSS - 2,000 bp` through `TSS + 500 bp`, clipped to chromosome bounds.
- Sensitivity windows: `-1,000/+100` and `-5,000/+1,000` bp.
- An overlap requires at least one shared base under half-open semantics. Boundary-touching intervals do not overlap.
- A region overlapping multiple promoters yields one row per region–gene pair. A gene linked to multiple regions retains every link; gene summaries aggregate only in separate derived tables.
- Overlapping genes are not collapsed. Ambiguity is represented by repeated links and `multi_gene_region=true`.

## Assay-to-promoter rules

- **CpGs:** require configurable minimum coverage, default 5 reads per sample and presence in at least the number of samples required by the contrast. Preserve methylated and total counts. Promoter methylation is coverage-weighted: `sum(methylated_count) / sum(total_count)`. Also report CpG count, covered bases, median CpG proportion, and coverage. Differential region statistics take precedence over post-hoc averaging when a validated DMR method was run.
- **ATAC:** a consensus peak is promoter-associated on at least 1 bp overlap. Counts come from the frozen consensus peak set, not sample-specific peak presence/absence. Peak-level rows are retained; promoter summaries report count of linked peaks and a declared weighted or strongest-effect summary.
- **H3K27ac:** the same overlap policy applies to a frozen consensus peak set. Input-control identity and provenance are mandatory for raw ChIP processing. Signal-only processed inputs may support descriptive integration but not count-based differential enrichment.

## Distal candidate links

Promoter overlap is the highest-confidence coordinate class in v0.1. Distal links are optional and must identify one method:

1. `nearest_tss_candidate`: nearest eligible TSS within a default maximum of 100 kb; ties produce all tied genes.
2. `distance_window_candidate`: every eligible TSS within a configurable distance, default 100 kb.
3. `external_map`: an externally supplied, versioned enhancer–gene map with its original score and provenance.

Nearest TSS is never labeled the true target. Regions outside the selected distance with no external link remain unlinked. Contact-informed and ABC-style links are deferred.

## Direction and aggregation

- RNA: positive effect means higher expression in numerator/test versus denominator/reference.
- ATAC: positive means greater accessibility.
- H3K27ac: positive means greater enrichment.
- Methylation: positive means greater methylation; report as proportion-point change in addition to model scale where applicable.
- `concordant_activation_association` and `concordant_repression_association` are descriptive composite labels, not mechanistic conclusions.
- A gene with mixed directions across its linked regions is `discordant_or_heterogeneous`; the workflow must not select only the most concordant region.

## Missing modalities

Integration is an outer join over valid gene–region links. Modality-specific values remain `NA`; Boolean flags such as `rna_missing` identify absence. Zero is reserved for a measured numeric zero. A composite evidence class states exactly which modalities contributed and whether each contribution is observed, inferred, predicted, or independently validated.

## Processed-input provenance

Every external processed input preserves, where available: `source_genome_build`, `source_annotation`, `source_processing_pipeline`, `source_file_accession`, `source_file_checksum`, `coordinate_normalization_status`, and `gene_id_mapping_status`. Original coordinates and identifiers must remain recoverable after normalization or harmonization.

`UNKNOWN` means the source did not provide a value. `UNRESOLVED` means available metadata are conflicting, ambiguous, or not yet verified. Neither value may be replaced by a production default. In particular, choosing GENCODE 50 for production does not relabel the source annotation of an external file.

## Integrated table

Required columns are:

`gene_id`, `gene_id_versioned`, `original_gene_id`, `gene_name`, `region_id`, `region_type`, `chromosome`, `start`, `end`, `strand`, `coordinate_system`, `genome_build`, `source_genome_build`, `source_annotation`, `source_processing_pipeline`, `source_file_accession`, `source_file_checksum`, `coordinate_normalization_status`, `gene_id_mapping_status`, `link_method`, `link_distance_bp`, `rna_effect_size`, `rna_p_value`, `rna_adjusted_p_value`, `atac_effect_size`, `atac_p_value`, `atac_adjusted_p_value`, `methylation_effect_size`, `methylation_p_value`, `methylation_adjusted_p_value`, `h3k27ac_effect_size`, `h3k27ac_p_value`, `h3k27ac_adjusted_p_value`, `evidence_class`, `rna_missing`, `atac_missing`, `methylation_missing`, `h3k27ac_missing`, `multi_gene_region`, `analysis_state`, `provenance_id`.

`provenance_id` resolves to input checksums, software/container identifiers, reference bundle, parameters, contrast, and transformation history. Sensitivity analyses create separately named result sets; they do not overwrite the primary result.
