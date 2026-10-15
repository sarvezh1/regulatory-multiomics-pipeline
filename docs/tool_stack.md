# v0.1 tool stack

The RNA stack is frozen from RNA validation execution. The ATAC stack is frozen from ATAC Docker/Linux-amd64 validation. Later-modality sections retain the initial audited selections. Exact immutable container references are in `assets/tool_containers.tsv`.

## Common and reporting

| Role | Selection | Validated pin | Rationale |
|---|---|---:|---|
| raw-read QC | FastQC | 0.12.1 | Stable assay-agnostic metrics and broad MultiQC support. |
| read preprocessing | fastp | 0.23.4 | Paired-read aware, JSON/HTML metrics, and stable container availability. |
| alignment utilities | samtools | 1.20 | Mature BAM operations in the executed branch. |
| interval operations | BEDTools | 2.31.1 | Transparent overlap operations; custom scripts still enforce coordinate contracts. |
| aggregate QC | MultiQC | 1.35 | Deterministic local parsing with remote upload disabled; publication uses a static summary and curated parsed TSV tables. |
| RNA statistical environment | R / Bioconductor | R 4.5.2 / Bioconductor 3.22 | Executed release family for DESeq2 1.50.2. |

Future modalities must select one coherent Bioconductor release family when implemented. They may not silently inherit a newer host installation or mix incompatible release families. ATAC deliberately remains on the RNA R 4.5.2/Bioconductor 3.22 family.

## RNA-seq

- FastQC 0.12.1, fastp 0.23.4, STAR 2.7.11b, samtools 1.20, featureCounts/Subread 2.1.1, DESeq2 1.50.2, and MultiQC 1.35.
- STAR provides splice-aware genomic alignment. featureCounts generates exon-grouped integer `gene_id` counts from the matching GTF.
- Raw counts feed DESeq2. Size-factor-normalized counts are descriptive only.

Salmon is not part of the RNA production branch because transcript-level inference is not required for promoter-oriented gene-level integration. It may be evaluated later without changing the canonical gene-count contract. edgeR remains an evaluated alternative, not a second default.

## ATAC-seq

- FastQC 0.12.1, fastp 0.23.4, Bowtie2 2.5.5, samtools 1.20, MACS3 3.0.5, ataqv 1.3.1, deepTools 3.5.6, BEDTools 2.31.1, and MultiQC 1.35.
- Statistical environment: R 4.5.2, Bioconductor 3.22, DiffBind 3.20.0, and DESeq2 1.50.2.
- Raw support is deliberately paired-end only. MACS3 uses observed fragments in BAMPE mode. Picard, motif analysis, footprinting, and IDR are not part of this branch.
- The canonical integer region matrix is the DA boundary. DiffBind is loaded and version-validated in the production statistical environment; the precomputed common-region matrix is fitted through its DESeq2 count-model backend without BAM recounting.

The direct matrix boundary provides identical inference behavior for raw-derived and validated externally processed counts. It does not claim that external counts were produced by this workflow.

## WGBS

- FastQC 0.12.1 and Trim Galore 2.3.0 perform raw QC and conservative paired trimming (`--quality 20 --length 20`); RRBS flags are not used.
- Bismark Rust suite 3.1.0 with its bundled Bowtie2 2.5.5 performs deterministic reference preparation, faithful/default short-read alignment, per-library paired deduplication, and methylation extraction. The compatibility methylation-extractor report identifies its classic extractor component as 0.25.1.
- R 4.5.2, Bioconductor 3.22, dmrseq 1.30.0, and bsseq 1.46.0 provide the sole production region-level differential method. DSS is not implemented.
- MultiQC 1.35 aggregates supported WGBS reports. Python 3.13.13 runs canonical contract/provenance code; the dmrseq image supplies Python 3.12.7 only for its local wrapper task.

WGBS retains methylated and unmethylated counts and never substitutes percentages for counts. Bismark indexes are distinct from ordinary Bowtie2 indexes. WGBS validation built only the execution-only chr22 index.

## H3K27ac ChIP-seq (planned; not implemented)

- fastp 1.3.7, Bowtie2 2.5.4, samtools 1.22.1, MACS3 3.0.5, deepTools 3.5.6.
- DiffBind plus DESeq2 from Bioconductor 3.23 for replicated enrichment testing.
- ChIPseeker and GenomicRanges for annotation/integration.

Every ChIP library must resolve to a compatible input control for raw processing. Control-free processed signal is explicitly descriptive. Broad-versus-narrow peak behavior will be fixture-tested; the primary H3K27ac peak policy uses MACS3 broad peaks plus reproducible signal tracks and retains the exact call parameters.

## Integration implementation (planned; not implemented)

GenomicRanges and rtracklayer define in-memory genomic operations and import/export. `data.table` is allowed for large deterministic joins; selected tidyverse packages may be used in reports but are not required by core contracts. Small custom scripts will validate schemas, normalize identifiers, perform deterministic joins, and emit provenance. Their algorithms must be unit-tested and must not replace established differential models.

## Maintenance, containers, and footprint

The selected tools have maintained upstream or Bioconductor release channels and are available through BioContainers/Bioconda or official Bioconductor images. Each implemented branch must:

1. select one image per coherent tool group;
2. record registry, tag, digest, architecture, licenses, and package inventory;
3. run a compact compatibility test before accepting the lock;
4. reject floating tags such as `latest`;
5. keep reference-data checksums separate from image provenance.

All genome indexes must be built from the exact pinned GRCh38 primary-assembly FASTA defined by the production reference policy. GENCODE 50 supplies the production annotation but does not, by itself, identify the FASTA representation. An external `hg38` label is insufficient evidence that an index or processed file matches the production contig set.

Expected compressed image cache for the eventual complete workflow is 30–100 GB. Full-GRCh38 STAR indexing needs roughly 30–40 GB RAM and disk; full-scale WGBS remains the dominant compute and temporary-storage path. No full production Bismark index has been built locally.

## Sources

- [Bioconductor release compatibility](https://bioconductor.org/about/release-announcements/)
- [DESeq2 1.50.2 in Bioconductor 3.22](https://bioconductor.org/packages/3.22/bioc/html/DESeq2.html)
- [fastp releases](https://github.com/OpenGene/fastp/releases)
- [MACS changelog](https://github.com/macs3-project/MACS/blob/main/ChangeLog)
- [ENCODE portal data model and API](https://www.encodeproject.org/help/rest-api/)
- [Bismark documentation](https://felixkrueger.github.io/Bismark/)
- [dmrseq 1.30.0 in Bioconductor 3.22](https://bioconductor.org/packages/3.22/bioc/html/dmrseq.html)
