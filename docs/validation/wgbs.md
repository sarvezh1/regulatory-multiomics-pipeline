# WGBS validation record

## Validation boundary

This record validates the conventional paired-end WGBS implementation, contracts, guards, and a bounded public execution. The compact public run uses only GRCh38.p14 chromosome 22 and is not a genome-wide benchmark, biological analysis, conversion-efficiency study, or quality pass/fail assessment. No production GRCh38 Bismark index was built.

## Synthetic raw execution

The deterministic fixture contains two technical libraries for one biological sample, 34 pairs per library, deliberate coordinate duplicates, known methylated and unmethylated CpGs, and reads from both complementary strands. The production path completed through FastQC, Trim Galore, Bismark reference preparation/alignment/deduplication/extraction, canonical CpG construction, aggregation, MultiQC, and provenance.

Both libraries mapped 34/34 pairs. Each reported duplicate removal, recovered reference CpGs, produced M-bias output, and satisfied exact count identity. Canonical sample loci showed `bottom,top` provenance and both contributing libraries. Processing, extraction, and count-contract statuses were `COMPLETE`; differential analysis alone was `NOT_RUN_DISABLED`.

An explicit zero-alignment fixture completed processing with zero unique alignments, `wgbs_methylation_extraction = COMPLETE_EMPTY_RESULT`, an empty valid count contract, and disabled differential status.

## Synthetic processed and differential execution

The replicated fixture has three controls and three treatments with unchanged, hypermethylated, hypomethylated, low-coverage, missing, isolated, and multi-CpG stretches. Ninety-eight CpGs passed the frozen rule. dmrseq 1.30.0 returned three complete candidate-region rows; the intended hypermethylated region had numerator delta +0.70 and the intended hypomethylated region -0.70. Missing coverage remained missing. No returned region had q-value below 0.05, demonstrating that a valid zero-significant-result analysis remains `COMPLETE`.

Guards produced the intended states:

- 1 versus 1: `NOT_RUN_INSUFFICIENT_REPLICATION`
- missing contrast: `NOT_RUN_NO_CONTRAST`
- rank-deficient design: `NOT_RUN_INVALID_DESIGN`
- percentage-only input with inference requested: `NOT_AVAILABLE_INPUT`
- all loci below coverage five: `COMPLETE_EMPTY_RESULT`
- zero covered loci: `COMPLETE_EMPTY_RESULT` with a distinct reason

Validation rejects missing R2, unknown raw directionality, unsupported protocols, mixed processed tiers, incompatible unresolved builds, negative counts, inconsistent coverage, fractions outside [0,1], invalid coordinates, and zero coverage presented as zero methylation. Reference-aware plus/minus merging and explicit non-reference-CpG exclusion have unit coverage.

## Genuine public execution

GSE152146 was not used for raw validation because the primary publication documents NEBNext Enzymatic Methyl-seq rather than conventional bisulfite sequencing. This protocol distinction was discovered during implementation and is preserved here without changing the frozen earlier audit.

The replacement is the already audited ENCODE K562 conventional WGBS experiment [ENCSR765JPC](https://www.encodeproject.org/experiments/ENCSR765JPC/): GSM2308597, SRX2157008, SRR4235744, SAMN05733812, biological replicate 1, paired Illumina HiSeq 2500. The [primary HAIB protocol](https://www.encodeproject.org/documents/6e461240-5af7-462c-843f-3dae804f957a/@@download/attachment/HAIB_Myers_WGBSprotocol_2014-10-30.pdf) ligates methylated paired-end adapters before bisulfite conversion and then amplifies with the corresponding paired-end primers; this supports the explicit directional classification. ENA reports source mate FASTQs of 66,943,138,587 and 64,011,469,703 bytes with MD5 values recorded in `tests/integration/genuine_wgbs_run_metadata.tsv`.

The pinned SRA Toolkit 3.4.1 container streamed only the first 250,000 complete spots from the public SRA object. Deterministically compressed local mates total 32,337,377 bytes. Execution used the existing chr22 FASTA (`sha256:dbfe87541379f14d76ba631000c76e1fdb80c5f2f3c4df14c84395ef5378dedc`) and matching GENCODE v50 chr22 GTF.

Diagnostics were: 250,000 input pairs; 189,997 pairs retained after trimming; 3,303 unique chr22 alignments (1.7% compact-reference mapping); one duplicate removed (0.03%); 5,721 canonical observed CpGs; 10 calls lacking a reference CpG at the reported coordinate explicitly excluded; 14 CpGs at coverage >=5; mean/median observed coverage 1.0736/1; observed CpG methylation 0.2976; CHG 0.0104; CHH 0.0097. M-bias was generated. Conversion-control efficiency was `NOT_AVAILABLE_INPUT` because the compact reference does not include the protocol's lambda control contig. These values are execution diagnostics only.

The first full execution took about 3 minutes 39 seconds including index construction. The compact-reference scope was then explicitly recorded as `compact_chr22_execution`, with annotation release `GENCODE 50`, rather than inheriting the synthetic test-profile label. A subsequent cache-resumed QC correction completed in about 30 seconds. Peak recorded RSS was approximately 992 MB for the disabled differential wrapper and 398 MB for compact Bismark reference preparation. Published genuine-run output occupied approximately 292 MB.

## Executed versions and containers

Validated versions are Nextflow 26.04.6; FastQC 0.12.1; Trim Galore 2.3.0; Bismark 3.1.0; Bowtie2 2.5.5; R 4.5.2; Bioconductor 3.22; dmrseq 1.30.0; bsseq 1.46.0; MultiQC 1.35; and contract Python 3.13.13. The WGBS path does not independently invoke samtools. Exact image references and resolved digests are in `assets/tool_containers.tsv`.

Processed real-data adapter validation was not added because the available public files were hundreds of megabytes and the synthetic count-preserving route already exercised the same adapter without unnecessary download.

## Regression and closure checks

The complete WGBS offline runner completed successfully, including raw, processed, differential, guard, empty-result, and zero-alignment scenarios; its 14 WGBS contract unit tests also passed. The established RNA and ATAC offline runners subsequently completed successfully with their output assertions and shared contract tests. The WGBS output assertion now also verifies the genuine fixture's retained-pair count and compact-reference provenance.

The retained WGBS fixtures and published test outputs occupy 386,972,689 bytes, including 32,337,596 bytes for the bounded genuine fixture and 291,667,065 bytes for its published output. The shared Nextflow work cache occupied approximately 4.05 GB at closure. D: retained 91,083,694,080 bytes free, comfortably above the 50 GB stop boundary. No full cohort or full-GRCh38 Bismark index was downloaded or built.
