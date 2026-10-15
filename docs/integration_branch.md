# Regulatory integration branch

## Scope

The integration branch combines standardized RNA-seq, ATAC-seq, conventional WGBS, and H3K27ac ChIP-seq contracts. It reports observed measurements, upstream assay-specific inference, deterministic coordinate relationships, and explicitly labelled candidate links. It does not infer mechanism or causality.

Two modes are distinct:

- `baseline` associates descriptive measurements and genomic features without differential labels.
- `differential` accepts only upstream results whose component status is `COMPLETE` and whose explicit variable, numerator, and denominator agree across modalities.
- `both` produces the descriptive layer and includes compatible differential evidence.

Partial-modality manifests are valid. Each modality is recorded as available or unavailable and as tested or not tested; absent or untested evidence is never converted to a zero effect.

## Invocation and manifest

```bash
nextflow run main.nf -profile docker \
  --assay INTEGRATION \
  --integration_manifest integration_manifest.tsv \
  --fasta GRCh38.primary_assembly.fa \
  --gtf gencode.v50.annotation.gtf \
  --integration_mode both \
  --integration_match_type EXACT_SAMPLE_MATCHED
```

The manifest schema is defined in `schemas/integration_manifest.tsv`. Data types are exact, not fuzzy. Coordinate-bearing inputs require an explicit compatible genome build and a resolved coordinate-normalization state. The accepted biological correspondence classes are `EXACT_SAMPLE_MATCHED`, `CONDITION_MATCHED`, `REFERENCE_COMPARABLE`, and `UNMATCHED`.

## Reference, identifiers, and promoters

The production coordinate space is the GRCh38 primary assembly on chromosomes 1–22, X, and Y with GENCODE release 50 annotation. Mitochondrial, alternate, patch, haplotype, unlocalized, and unplaced sequence is excluded from default integration.

Cross-modality gene joins use the base Ensembl gene identifier produced by deterministically removing only a terminal numeric version. Original identifiers and mapping records remain available. Symbols are descriptive and are not primary keys.

The integration TSS resource retains every distinct GENCODE transcript TSS and collapses only exact duplicate coordinates for the same gene and strand. Primary promoters are strand-aware TSS -2000/+500 intervals, represented as zero-based half-open coordinates and clipped at chromosome boundaries. Different TSSs remain different promoters, including overlapping and alternative promoters.

## Overlap and aggregation rules

Genomic overlap requires a positive base-pair intersection. Book-ended intervals do not overlap. Detailed ATAC–promoter, H3K27ac–promoter, CpG–promoter, and DMR–promoter tables preserve one-to-many relationships.

Promoter methylation is coverage weighted. For count-preserving inputs it is `sum(methylated) / sum(coverage)`; uncovered loci remain missing. Percentage-plus-coverage public inputs are explicitly marked descriptive, retain unknown methylated and unmethylated read counts, and cannot support differential inference.

Promoter summaries report counts of associated regions and directions while retaining source identifiers in the detailed tables. Gene summaries retain promoter counts and report `PROMOTER_HETEROGENEITY` when promoter patterns disagree.

## Differential evidence and patterns

The integration layer does not rerun assay-level differential analysis. Positive effects retain upstream meanings: greater numerator expression, accessibility, methylation, or H3K27ac enrichment. The default evidence threshold is adjusted p-value or q-value at or below 0.05 and is configurable before analysis.

Literal states are `UP`, `DOWN`, `HYPERMETHYLATED`, `HYPOMETHYLATED`, `NO_SIGNIFICANT_CHANGE`, `NOT_TESTED`, and `NOT_AVAILABLE`. Pattern strings concatenate these states without mechanistic labels. Discordant patterns are retained as results.

No cross-assay p-value combination or aggregate regulatory score is calculated.

## Co-localization and distal candidates

Direct ATAC/H3K27ac overlap is labelled `ATAC_H3K27AC_COLOCATED`. This is biochemical co-localization, not proof of enhancer function.

Distal candidate regions are anchored on ATAC when ATAC is present; otherwise H3K27ac-only anchors are permitted. The primary rule links a non-promoter region to the nearest annotated TSS within 100 kb. Exact distance ties retain every tied TSS. Links are labelled `NEAREST_TSS_CANDIDATE` and do not imply contact or a true target.

Optional external enhancer–gene links retain their own source, version, genome build, evidence type, and validation type. External predicted and experimentally validated relationships remain distinct from nearest-TSS candidates. ABC and contact-map computation are deferred because linear proximity alone does not establish physical enhancer–promoter contact.

## Sensitivity, output, and statuses

The primary promoter window (-2000/+500) is compared with predefined -1000/+500 and -2000/+1000 windows. The 100 kb distal limit is compared with 50 kb and 250 kb. The workflow reports association, candidate-link, and pattern changes without selecting the setting that maximizes concordance.

Outputs under `integration/` include reference tables, compatibility reports, detailed relationships, promoter and gene evidence tables, sensitivity summaries, plot source data and SVG, component statuses, and provenance. Status components are `integration_reference`, `integration_promoter`, `integration_differential`, and `integration_distal` and use the shared project vocabulary.

Descriptive correlations are not implemented. Such correlations would require demonstrated exact sample matching and adequate matched sample size.
