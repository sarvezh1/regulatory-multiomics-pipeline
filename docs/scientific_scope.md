# v0.1 scientific scope

Author: Sarvesh Galgale
Status: frozen for implementation review
Specification date: 2026-10-12

## Purpose and claim boundary

The workflow will integrate human bulk short-read RNA sequencing, ATAC-seq, conventional WGBS, and H3K27ac profiling to examine associations among transcription, accessibility, DNA methylation, and H3K27ac-associated regulatory activity. It will keep four evidence levels distinct:

1. **Observed measurement**: counts, coverage, methylated counts, signal, or called regions directly supported by an assay.
2. **Statistically inferred change**: an effect estimate and uncertainty from a valid replicated contrast.
3. **Predicted regulatory association**: a coordinate- or map-based region-to-gene link.
4. **Experimentally validated relationship**: evidence explicitly supplied from an independent perturbation or targeted experiment.

Correlation, coordinate overlap, concordant direction, and nearest-gene assignment do not establish causation. Reports must use *associated*, *candidate*, or *predicted* unless independent validation is supplied.

## Supported in v0.1

- Homo sapiens bulk samples in the frozen production reference space defined below.
- Short-read RNA-seq, ATAC-seq, conventional WGBS, and H3K27ac ChIP-seq.
- Raw-input and validated processed-input entry points.
- Each modality alone, any valid subset, and complete four-modality integration.
- Assay-specific QC and, only where design is valid, assay-specific differential analysis.
- Promoter integration and cautious distal candidate linking.
- Deterministic tables, reports, software provenance, reference checksums, and explicit analysis states.
- Sensitivity analysis for promoter windows, peak sets, CpG coverage, and link rules.

## Frozen production reference policy

- Organism: Homo sapiens.
- Annotation: GENCODE release 50. GENCODE v50 is based on GRCh38.p14.
- Coordinate representation: the GRCh38 primary-assembly FASTA selected for production, with its exact filename, source URL, checksum, sequence dictionary, and chromosome sizes recorded.
- Default cross-modality integration space: primary chromosomes `chr1`â€“`chr22`, `chrX`, and `chrY`. Mitochondrial reads may be retained for assay QC but `chrM` is excluded from default nuclear regulatory integration.
- Assembly patches, alternate loci, haplotypes, unlocalized sequences, and unplaced sequences are excluded from default cross-modality integration. Supporting any of them requires an explicit future specification and compatible references for every modality.

GENCODE 50 is an **annotation release**; GRCh38.p14 identifies the genome assembly patch release on which it was annotated. Neither term alone identifies the exact FASTA or coordinate representation. The production bundle therefore pins the chosen primary-assembly FASTA independently from the GENCODE 50 GTF. A source dataset labelled only `hg38` is not assumed to use the full GRCh38.p14 all-regions assembly, the selected production FASTA, or GENCODE 50 annotation.

Inputs declaring hg19/GRCh37 are rejected unless a separately validated coordinate-normalization step is explicitly requested. Processed `hg38`/GRCh38 inputs require sequence-name and contig compatibility checks before integration. Silent lift-over, silent contig dropping, and mixed coordinate spaces are forbidden.

## Out of scope in v0.1

- Hi-C processing or contact-informed/ABC linking.
- Single-cell or spatial assays; proteomics.
- Direct RNA and required long-read RNA input.
- Non-human and de novo/no-reference analysis.
- CUT&Tag production support.
- Multi-mark ChromHMM state reconstruction.
- Causal claims and black-box machine-learning classification.

## H3K27ac decision

**Decision A: ChIP-seq is the only production/default H3K27ac branch in v0.1; CUT&Tag is deferred.** Both anchor datasets use H3K27ac ChIP-seq. ChIP-seq requires explicit input-control provenance, while CUT&Tag has different background, duplicate, fragment, and peak-calling behavior. Supporting both before a CUT&Tag validation dataset exists would double assay-policy surface without validating it. A later version may add CUT&Tag through the same output contract after compact and biological validation.

## Frozen biological questions

The workflow may answer:

- Q1: Which genes are differentially expressed?
- Q2: Which genomic regions change accessibility?
- Q3: Which CpGs or regions change DNA methylation?
- Q4: Which regions change H3K27ac enrichment?
- Q5: Which genes have promoter-associated regulatory changes concordant with transcription?
- Q6: Which genes show discordant regulatory patterns?
- Q7: Which distal candidate regulatory elements can be linked cautiously to genes?
- Q8: How sensitive are conclusions to promoter definition, peak set, CpG coverage threshold, and regulatory-region linking rule?

Patterns such as hypermethylation with decreased accessibility, H3K27ac, and expressionâ€”or the reverseâ€”are reportable observations, not universal laws. A question is answered only for modalities with valid inputs and designs.

## Statistical design requirements

- Biological and technical replicates remain distinct. Technical replicates may be merged only within the same biological library/sample under a recorded rule; they never increase biological `n`.
- Every differential run requires a named contrast, recorded model formula/design matrix, estimable coefficients, and independent biological replication in each contrasted group. Default minimum: at least two biological replicates per group; three or more is strongly preferred.
- Batch or confounders are included only when scientifically justified and the design remains full rank. Paired subjects use an explicit subject blocking term.
- Effect sizes, standard errors where available, raw p-values, adjusted p-values, tested-universe size, and multiple-testing method are emitted.
- An unreplicated assay may produce descriptive signal and coordinate integration but no inferential differential result.
- Missing and untestable values remain `NA`; they are never replaced with zero.

Every analysis unit ends in exactly one state:

- `COMPLETE`
- `COMPLETE_EMPTY_RESULT`
- `NOT_RUN_DISABLED`
- `NOT_RUN_NO_CONTRAST`
- `NOT_RUN_INSUFFICIENT_REPLICATION`
- `NOT_RUN_INVALID_DESIGN`
- `NOT_AVAILABLE_INPUT`
- `FAILED`

`COMPLETE_EMPTY_RESULT` means a valid test completed and returned no features meeting a declared filter. It is not interchangeable with any `NOT_RUN_*` state.

## Sources

- [GENCODE human release history](https://www.gencodegenes.org/human/releases.html)
- [ENCODE K562 reference epigenome ENCSR612NLL](https://www.encodeproject.org/reference-epigenomes/ENCSR612NLL/)
- [GEO GSE152146](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152146)
