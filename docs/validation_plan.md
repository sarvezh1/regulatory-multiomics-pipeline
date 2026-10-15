# Validation and execution plan

## Four validation layers

### 1. Structural and synthetic

Tiny deterministic fixtures will test samplesheet types and uniqueness; paired-read rules; malformed metadata rejection; genome-build mismatch; interval boundaries; 1-based to 0-based conversion; chromosome clipping; overlapping genes; multiple TSSs; promoter and distal links; missing modalities; `NA` preservation; effect direction; evidence classification; technical-versus-biological replication; invalid contrasts; and every analysis-state value.

Expected values are hand-auditable. Tests establish software behavior, not biological truth.

### 2. Compact genuine data

Each raw branch must eventually execute on small genuine or reproducibly downsampled reads: RNA, ATAC, WGBS, and H3K27ac ChIP with input control. Fixtures must retain real read-pair structure and assay characteristics. Downsampling commands, source accessions, and checksums are recorded. These tests prove execution and file contracts, not sensitivity.

### 3. Processed-data integration

Use released processed subsets from ENCSR612NLL for cross-assay reference, coordinate, and contract integration. These are independent K562 experiments, not matched aliquots. Use GSE152146 expression tables, ATAC/H3K27ac bigWigs, and Bismark coverage selectively; do not fetch its 5.2 GB aggregate archive for routine validation. Validate source-provenance preservation, exact reference compatibility, coordinate normalization, identifier mapping, missingness, and deterministic joins.

### 4. Biological validation

Use GSE152146 only within its design limits:

- four-modality baseline/integration role: untreated RNA + ATAC + WGBS + H3K27ac, without treating the assays as replicate-matched;
- partial-modality differential role: only modalities and contrasts with legitimate independent biological replication may undergo differential inference;
- replicated RNA+ATAC: untreated versus 12 h TNF, and untreated versus recovery, are legitimate two-modality contrasts;
- H3K27ac in the same states is descriptive because each state has one library;
- WGBS supports descriptive untreated and 12-day comparisons but no replicated differential inference;
- no nontrivial four-modality differential contrast exists.

GSE152146 must never be presented as a replicated four-modality differential benchmark.

Agreement with the publication is a review target, not a unit-test oracle. Unexpected valid results must remain possible.

## Execution tiers and storage

| Tier | Purpose | Expected local footprint |
|---|---|---:|
| 1 | synthetic logic and CI | <2 GB including outputs |
| 2 | compact genuine raw paths | 20–80 GB transient, <15 GB retained |
| 3 | public processed integration | 15–100 GB selected; full raw cohorts on external storage/HPC |

Planning ranges for a modest two-condition, three-replicate design are intentionally conservative:

| Component | Compressed input | Peak working/outputs |
|---|---:|---:|
| RNA-seq | 30–90 GB | 100–300 GB |
| ATAC-seq | 20–80 GB | 80–250 GB |
| H3K27ac plus inputs | 20–100 GB | 100–350 GB |
| WGBS | 300 GB–1.5 TB | 1–5 TB |
| GRCh38 indexes | — | STAR 30–40 GB; Bowtie2/Bismark 15–30 GB; shared references 10–25 GB |
| Container cache | — | 30–100 GB |

Nextflow-style work directories can reach roughly 2–4 times retained deliverables because staged inputs and intermediates coexist. Actual size depends on depth, replicate count, compression, and retention; a preflight estimator must use declared file sizes.

Future execution behavior must provide independent modality branches, processed-data entry points, resume, publish by link where safe, no unnecessary copies, optional post-success cleanup, configurable intermediate retention, scratch-directory selection, and disk-space preflight. WGBS Tier 3 is not assumed feasible on the current workstation.

## Acceptance gates for RNA implementation

- Frozen schemas validate all anchor metadata examples.
- Reference bundle and container digests are recorded before execution.
- The production FASTA representation is pinned separately from GENCODE 50, and its default integration contigs are verified as primary chromosomes only.
- Every external processed input records source build, annotation, pipeline, accession, checksum, coordinate-normalization status, and gene-ID mapping status using explicit `UNKNOWN`/`UNRESOLVED` values where necessary.
- Exact H3K27ac input-control identities and required ENCODE raw-file metadata are resolved before raw H3K27ac validation. This gate is scoped to raw H3K27ac and does not block unrelated modalities.
- Each process has a declared output contract and terminal state.
- Tests explicitly cover insufficient replication and invalid designs.
- No large public dataset is required to run routine tests.

## RNA execution record

The RNA-specific structural, synthetic, guard, compact genuine-data, provenance, and reporting gates were executed. Exact fixtures, checksums, measured outputs, versions, and limitations are recorded in `docs/validation/rna.md`. Later-modality tests in this plan remain future work and were not started.
