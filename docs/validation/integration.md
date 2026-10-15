# Integration validation

## Deterministic validation

The synthetic fixture exercises all four modalities plus partial and contrast-mismatch runs. It includes concordant and discordant directions, untested versus nonsignificant evidence, missing modalities, alternative promoters, multiple regions and DMRs per promoter, shared promoters, chromosome-boundary clipping, versioned gene identifiers, ATAC/H3K27ac co-localization, tied nearest-TSS candidates, an out-of-range distal region, and an externally validated relationship.

Assertions cover identifier preservation, strand-aware coordinates, half-open boundaries, detailed links, count-weighted methylation, pattern strings, contrast rejection, partial integration, tie preservation, distance limits, sensitivity outputs, provenance, and component-scoped statuses.

## Genuine descriptive validation

The genuine validation uses independent released K562 ENCODE experiments and is classified `REFERENCE_COMPARABLE`, not sample matched:

| Modality | Experiment | File | Use |
|---|---|---|---|
| RNA | ENCSR109IQO | ENCFF611MXW | GENCODE V29 gene quantification; ENSG rows and descriptive TPM only |
| ATAC | ENCSR868FGK | ENCFF057UYP | GRCh38 replicated narrowPeak, restricted to chr22 |
| WGBS | ENCSR765JPC | ENCFF750NZN | GRCh38 CpG bigBed, remotely range-read for chr22 |
| H3K27ac | ENCSR000AKP | ENCFF038DDS | GRCh38 replicated narrowPeak, restricted to chr22 |

The WGBS bigBed exposes read coverage and rounded methylation percentage rather than exact methylated/unmethylated counts. Its adapter therefore produces descriptive CpG evidence only, preserves count fields as `NA`, and does not permit differential inference. The full 1.14 GB file was not downloaded.

The execution-only chr22 reference is GRCh38.p14 primary assembly with the matching GENCODE v50 chr22 annotation. The genuine run validates coordinate compatibility, promoter annotation, four-modality contract interoperability, interval relationships, partial quantitative availability, and missing-data semantics. It does not validate differential regulation, within-sample correlation, causality, enhancer targets, or biological mechanisms.

The validated run produced 8,676 promoters for 1,747 genes, 7,133 ATAC–promoter links, 6,457 H3K27ac–promoter links, 694,225 CpG–promoter links, 1,050 ATAC/H3K27ac co-localizations, and 1,884 nearest-TSS candidate links. `integration_differential` was `NOT_RUN_DISABLED`. The integration task completed in 94.8 seconds with 992,272,384 bytes peak RSS; reference derivation completed in 63 seconds with 1,103,396,864 bytes peak RSS.

Optional genuine GSE152146 RNA/ATAC differential integration is deferred. The synthetic fixture fully validates contrast-aware partial differential behavior without acquiring additional processed data; unreplicated GSE152146 WGBS and H3K27ac remain excluded from inference.

## Tools

Reference derivation used R 4.5.2, Bioconductor 3.22, GenomicRanges 1.62.1, IRanges 2.44.0, S4Vectors 0.48.0, and rtracklayer 1.70.1. Contract integration used Python 3.13.13. Nextflow 26.04.6 orchestrated the validated Docker Linux/amd64 execution. Exact immutable image references are recorded in `assets/tool_containers.tsv`.

RNA, ATAC, WGBS, and H3K27ac offline regression suites passed after the shared entry point and configuration changes. Integration fixtures and retained published validation outputs occupy 564,671,325 bytes (approximately 565 MB decimal); D: retained approximately 89.82 GB free at closure.

## Limitations

- Baseline abundance and signal values from independent projects are descriptive and are not directly compared across datasets.
- Nearest-TSS distal links are predictions based on linear distance.
- No contact-informed linking, causal inference, cross-assay omnibus test, or opaque combined score is performed.
- Public processed inputs may have older annotations; source annotation is retained rather than relabelled as GENCODE 50.
