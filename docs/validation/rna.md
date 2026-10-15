# RNA validation record

Author: Sarvesh Galgale

## Deterministic validation

The synthetic raw fixture contains one single-end library, a 3 kb `chr22` reference, two annotated genes, and 24 uniquely placeable 50-base reads. The completed workflow produced a BAM/BAI, 24 mapped primary reads, and exact counts of 12 for each expected gene. Component-scoped status records reported RNA processing, counting, and QC as `COMPLETE`; only the explicitly disabled differential component reported `NOT_RUN_DISABLED`.

The processed fixture contains three control and three treatment biological samples, 100 genes, deterministic heterogeneity, strong positive and negative effects, unchanged genes, low counts, and an all-zero gene. DESeq2 completed with the explicit contrast `condition,treatment,control`; intended effects had the correct sign, all 100 genes remained in the complete table, and the zero-count gene retained `NA` statistics.

Guard runs returned:

| Case | Expected and observed state |
|---|---|
| 1 control versus 1 treatment | `NOT_RUN_INSUFFICIENT_REPLICATION` |
| inference enabled without contrast | `NOT_RUN_NO_CONTRAST` |
| `~ batch + condition` with one batch level | `NOT_RUN_INVALID_DESIGN` |
| malformed paired metadata | validation failure |
| simultaneous raw and processed input | validation failure |

The technical-replicate fixture has seven libraries but six biological samples. Two libraries for `control_1` were summed before inference; the manifest records both, the combined matrix contains one `control_1` column, and DESeq2 reports independent biological n=6.

## Compact genuine execution

Source run: GSE152146 / PRJNA638409, GSM4715082, SAMN15760175, SRX8910158, SRR12413947. ENA records 65,070,148 source reads, BGISEQ, RNA-Seq, transcriptomic source, and paired library layout. Only the first 100,000 SRA spots were extracted with `fastq-dump -X 100000` through the official range-enabled NCBI object. The SRA record contains a zero-length technical mate and one usable 100-base biological read per spot, so the bounded fixture is represented truthfully as a single-read workflow input. The retained gzip is 9,610,295 bytes with SHA-256 `71d4da172b7818babb435830912fecf2677c11cbf3699d56a597cc16169263d4`.

The execution reference contains only GRCh38 chromosome 22 (`NC_000022.11`, renamed deterministically to `chr22`) and the matching GENCODE v50 `chr22` GTF subset. It is not a biological benchmark. Checksums are in `tests/integration/reference/checksums.tsv`; preparation commands are in `tests/integration/prepare_chr22_reference.sh`.

Observed execution diagnostics:

- 100,000 input records; 99,956 retained by fastp.
- 5,358 uniquely mapped reads, 908 multimappers, and 36 reads mapping to too many loci.
- 1,747 annotated genes in the count matrix; 356 had nonzero counts.
- 8 minutes 1 second for the first complete run; peak observed RSS was approximately 1.3 GB during STAR index generation.
- Retained genuine output: 457,455,456 bytes, including a 434,543,388-byte published compact STAR index.
- Entire local Nextflow work directory after all RNA runs and cache-invalidating verification reruns: about 877 MB; all test outputs: about 465 MB; compact reference/read fixtures: about 169 MB.

The 5.36% unique mapping rate is expected to be depressed by a chr22-only reference and bounded prefix. It is recorded only as an execution diagnostic and is not a production mapping benchmark or biological result.

## Validated runtime

| Component | Executed version |
|---|---:|
| Nextflow | 26.04.6 build 12646 |
| Java used by Nextflow | OpenJDK 17.0.20.1 |
| FastQC | 0.12.1 |
| fastp | 0.23.4 |
| STAR | 2.7.11b |
| samtools | 1.20 |
| featureCounts / Subread | 2.1.1 |
| R | 4.5.2 |
| Bioconductor release family | 3.22 |
| DESeq2 | 1.50.2 |
| MultiQC | 1.35 |
| contract runtime | Python 3.13.13 |

Container references and immutable registry digests are in `assets/tool_containers.tsv`. Bioconductor 3.22 and DESeq2 1.50.2 compatibility is recorded by the [Bioconductor 3.22 package page](https://bioconductor.org/packages/3.22/bioc/html/DESeq2.html). Docker on Linux/amd64 is the only container backend validated for RNA.

## Known limitations

- Full production GRCh38 indexing and publication-scale input were deliberately not tested locally.
- The genuine fixture validates single-file execution because of its SRA read structure; paired-end routing is schema-validated but should also receive a compact genuine conventional read-pair fixture in a later maintenance pass.
- Processed input accepts integer gene counts only; uncontrolled symbol mapping and non-integer abundance matrices are rejected.
- n=2 per group is permitted but weak; n>=3 remains preferred.
