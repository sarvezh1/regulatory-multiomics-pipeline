# ATAC validation record

Author: Sarvesh Galgale

## Foundation preflight

The raw RNA fixture already used component-scoped status: RNA processing/counting/QC were `COMPLETE`, while only differential expression was `NOT_RUN_DISABLED`. The documentation wording was clarified without redesigning the status model. The test entry point was hardened to fail on nonzero external-command exits. After the ATAC implementation and final reporting change, all six RNA workflow scenarios completed successfully, including the raw route, processed/DE route, technical-replicate case, and three inference guards. The outer regression wrapper expired while the final technical-replicate DESeq2 process was still running; that workflow subsequently completed with zero failed processes, a cache-confirmation rerun exited zero, output assertions passed, and all 15 shared Python contract tests passed.

## Implemented architecture

`main.nf` accepts only the explicit `RNA` or `ATAC` branch. Unsupported modality rows are rejected. The ATAC subworkflow has strict validation, reference/index/TSS preparation, library-level QC and alignment, deterministic fragment filtering, biological-sample merging, sample-level MACS3 peaks, CPM fragment-coverage bigWigs, ATAC QC, consensus construction, integer fragment counts, processed-matrix adaptation, differential guards/modeling, MultiQC, statuses, and provenance.

Raw ATAC is paired-end only. `sample_id` is biological material and `library_id` is a technical library. Technical libraries receive independent FastQC/fastp/alignment/filter metrics and combine within `sample_id` before peak calling/counting; they never add biological n.

## Frozen parameters and contracts

- Bowtie2: `--very-sensitive -X 2000 --no-mixed --no-discordant`.
- Filtering: name-sort, `samtools fixmate -m`, coordinate-sort, `samtools markdup -s`, primary properly paired reads, MAPQ >=30, mitochondrial exclusion, duplicate exclusion, then fragment-level GRCh38 blacklist exclusion.
- Duplicate marking occurs per technical library. Marked metrics and final nonduplicate BAMs are retained.
- Canonical BAM coordinates are unchanged. Derived insertions use +4/-5 only for TSS QC.
- MACS3: `-f BAMPE -q 0.01 --keep-dup all --call-summits`, narrow peaks, observed fragments, no single-end shift, no input control.
- Production effective genome size default: 2,913,022,398. Compact chr22 execution used 39,159,777 non-N bases, derived with `grep -v '^>' FASTA | tr -cd 'ACGTacgt' | wc -c`.
- Consensus regions merge overlap-connected peaks; book-ended intervals remain separate. IDs are `ATAC_` plus 16 hexadecimal SHA-256 characters from `chromosome:start-end`; contributing peaks remain traceable.
- Counting uses one paired fragment once per overlapping consensus region. Count columns are biological samples.
- Peak-only processed input is deferred. The implemented processed route requires coordinates, stable/source region IDs, biological-sample metadata, nonnegative integer counts, and source provenance.

The TSS resource uses unique deduplicated transcript TSS positions with strand and gene/transcript provenance in zero-based half-open BED. The compact GENCODE 50 chr22 GTF generated 8,664 unique positions. This QC resource is not the regulatory promoter window.

## Statistical route and states

DiffBind 3.20.0 constructs and validates the biological-sample/peak framework for the common region universe. Its DESeq2 1.50.2 count-model backend consumes the canonical integer matrix directly so raw-derived and processed routes have the same statistical boundary and no BAM recount is required. Supported safe formulas are `~ condition` and estimable additive covariate designs such as `~ batch + condition`. Contrasts use `variable,numerator,denominator`; positive log2 fold change means greater numerator accessibility.

Inference requires an explicit opt-in, explicit contrast, at least two independent biological samples per group, and a full-rank design. Complete tested-region tables preserve `NA`. The component states are `atac_processing`, `atac_peak_calling`, `atac_count_contract`, and `atac_differential`. A genuine zero-peak result is `COMPLETE_EMPTY_RESULT`, while a process crash remains a workflow failure.

## Synthetic validation

The paired-end raw fixture contains 950 input pairs, background fragments, two enriched loci, short and longer fragments, 50 mitochondrial fragments, deliberate duplicates, and blacklist overlaps. The complete real-tool path passed. It retained 950 mapped/proper fragments before filtering, recorded 39 duplicate fragments and 61 blacklist fragments, yielded 800 usable nuclear fragments, recovered two intended peaks/consensus regions with counts 352 and 351, and produced TSS, FRiP, fragment-length, complexity, ataqv, bigWig, MultiQC, provenance, and component-status outputs. The observed synthetic FRiP was 0.87875; it is a fixture assertion, not a universal threshold.

The processed DA fixture has three controls and three treatments across 100 regions. It preserves unchanged, positive, negative, low-count, and all-zero cases. DiffBind preserved all six biological samples; the DESeq2 backend completed, returned the intended effect directions, retained the complete 100-region table, and preserved unavailable all-zero statistics as `NA`.

Guard results were:

| Case | Observed result |
|---|---|
| 1 versus 1 | `NOT_RUN_INSUFFICIENT_REPLICATION` |
| missing contrast | `NOT_RUN_NO_CONTRAST` |
| rank-deficient design | `NOT_RUN_INVALID_DESIGN` |
| duplicate technical-library rows as processed samples | rejected |
| missing R2 | validation failure |
| raw plus processed conflict | validation failure |
| invalid assay | validation failure |
| malformed coordinates | validation failure |
| incompatible unnormalized build | validation failure |
| valid zero peaks | `COMPLETE_EMPTY_RESULT` |

The final offline suite passed 15 unit/contract tests, five ATAC Nextflow scenarios, and output assertions in 876.3 seconds. Genuine-data execution is excluded from the routine offline entry point. MultiQC 1.35 performs local parsing and aggregation with remote upload disabled. Retained QC publication uses a deterministic static HTML summary and a curated set of parsed TSV tables; temporary native report data are not published.

## Compact genuine execution

The source is GSE152146 / PRJNA638409, untreated ATAC biological replicates GSM4904180/SRX9486247/SRR13036031 and GSM4904181/SRX9486248/SRR13036032. Both are released paired-end Illumina NovaSeq 6000 runs. Full-run metadata, object sizes, ENA file sizes/checksums, and source URLs are in `tests/integration/genuine_atac_run_metadata.tsv`.

The bounded fixture was produced with SRA Tools 3.4.1 `fastq-dump <official SRA object URL> --split-files -X 250000`, followed by deterministic `gzip -n`. The initial 100,000-pair attempt yielded no chr22 peaks, so the same samples were increased within the predeclared 250,000-pair ceiling. Final FASTQ checksums are in `tests/integration/genuine_atac/checksums.tsv`; the four files total 70,095,347 bytes.

Execution used NC_000022.11 renamed to `chr22`, the matching GENCODE 50 chr22 GTF subset, a chr22 subset of ENCODE GRCh38 blacklist ENCFF356LFX, a deterministically built Bowtie2 index, and the derived TSS resource. Source/derived checksums are in `tests/integration/reference/checksums.tsv` and `atac_derived_checksums.tsv`; commands are in `prepare_chr22_reference.sh`. The compact reference is execution-only.

Observed diagnostics:

| Metric | SRR13036031 / GSM4904180 | SRR13036032 / GSM4904181 |
|---|---:|---:|
| input pairs | 250,000 | 250,000 |
| retained pairs | 247,604 | 247,968 |
| chr22 mapped/proper fragments | 11,673 | 11,138 |
| MAPQ >=30 fragments before exclusions | 2,755 | 2,648 |
| duplicate fragments | 867 | 912 |
| blacklist fragments | 21 | 29 |
| usable fragments | 2,447 | 2,342 |
| TSS enrichment | 3.3702 | 4.5523 |
| MACS3 q<=0.01 peaks | 0 | 0 |
| FRiP | 0.0 | 0.0 |

The chr22 reference intentionally omits `chrM`; mitochondrial availability is therefore `false` and mitochondrial fraction is `NA`, not a fabricated zero. MACS3 completed successfully with no q<=0.01 peaks, so peak calling is `COMPLETE_EMPTY_RESULT`; processing and the empty but valid count contract are `COMPLETE`; DA is `NOT_RUN_DISABLED`. These diagnostics are not genome-wide QC, treatment biology, mapping benchmarks, or biological pass/fail claims.

The 250,000-pair execution took approximately 12 minutes 37 seconds before cache-assisted closure reruns. Peak observed resident memory was about 545 MB (FastQC); Bowtie2 index construction used about 242 MB. The retained genuine ATAC output is 142,771,083 bytes. At final closure, the work directory is 3,184,050,578 bytes, all test outputs are 595,942,766 bytes, and all integration fixtures are 246,996,787 bytes. D: has 97,692,708,864 bytes free, versus approximately 100.1 GB at start. The incremental footprint remains far below the 15 GB preference and the 50 GB stop threshold.

## Validated runtime

| Component | Executed version |
|---|---:|
| Nextflow | 26.04.6 build 12646 |
| Java used by Nextflow | OpenJDK 17.0.20.1 |
| FastQC | 0.12.1 |
| fastp | 0.23.4 |
| Bowtie2 | 2.5.5 |
| samtools | 1.20 |
| MACS3 | 3.0.5 |
| ataqv | 1.3.1 |
| deepTools | 3.5.6 |
| BEDTools | 2.31.1 |
| R | 4.5.2 |
| Bioconductor release family | 3.22 |
| DiffBind | 3.20.0 |
| DESeq2 | 1.50.2 |
| MultiQC | 1.35 |

All scientific images use pinned tags and immutable digests recorded in `assets/tool_containers.tsv`. Docker on Linux/amd64 is the only validated backend.

## Implementation corrections and limitations

Corrections made during validation included hardened external-command exit handling; overlap-only versus book-ended consensus behavior; multi-sample version-file staging; empty-matrix QC and differential outputs; unavailable mitochondrial fractions on chr22; explicit effective-genome-size provenance; `NA` Jaccard for two empty peak sets; static allowlisted QC-report publication; and uncached lightweight provenance builders that re-hash files referenced indirectly by samplesheets without invalidating scientific-process resume. An initial genuine run also captured an interface change mid-execution; it was resumed and fully rerun after the interface was frozen.

Production-scale GRCh38 indexing and full-run ATAC processing were not performed. The genuine compact fixture has zero peaks because the bounded chr22-only data are sparse; nonempty real peaks are validated by the deterministic synthetic raw path, not claimed from GSE152146. Raw single-end ATAC, peak-only quantitative input, IDR, motifs, footprinting, and cross-modality integration remain out of scope.
