# H3K27ac validation

Author: Sarvesh Galgale

## Validation layers

Offline validation covers exact assay/target spelling, input-control relationships, single/paired layout guards, raw/processed conflicts, build compatibility, deterministic consensus IDs, book-ended intervals, empty-peak status, raw single-end and paired-end execution, processed 3-versus-3 differential enrichment, NA preservation, and inferential guards. Synthetic raw fixtures contain enrichment, background, uniquely named duplicates, chrM reads, and blacklist overlaps. The normal single-end MACS3 model and paired-end BAMPE path both complete.

The processed fixture returns the expected positive gain and negative loss, retains the all-zero region with NA inferential values, and reports `h3k27ac_differential = COMPLETE`. Guard executions return `NOT_RUN_INSUFFICIENT_REPLICATION`, `NOT_RUN_NO_CONTRAST`, and `NOT_RUN_INVALID_DESIGN` for 1-versus-1, absent contrast, and rank-deficient design respectively.

## Genuine compact execution

The technical anchor is released K562 H3K27ac experiment ENCSR000AKP. Its current ENCODE record identifies possible control ENCSR000AKY, a released K562 input-library experiment. The legacy primary metadata provides the explicit mapping `wgEncodeEH000043` to standard control `wgEncodeEH000052` and maps those records to GSM733656/SRX067407/SRR227385 and GSM733780/SRX067531/SRR227650. Both runs are single-ended 36-nt Illumina Genome Analyzer data. The current ENCODE experiment pages report one error and nine warnings for treatment and one error and four warnings for control. Exact experiment/file audit fields are recorded in `metadata/encode_h3k27ac_control_audit.tsv`; fields inaccessible from the current portal response are marked unresolved rather than inferred.

Exactly 500,000 spots per run were extracted with SRA Toolkit 3.4.1 from the public SRA objects and gzip-compressed deterministically. The retained SHA-256 values are recorded in `tests/integration/genuine_chip/checksums.tsv`. No complete experiment was downloaded.

Execution used the GRCh38.p14 primary-assembly chromosome 22 FASTA, matching GENCODE 50 chromosome 22 GTF, and the validated chr22 subset of ENCODE GRCh38 blacklist ENCFF356LFX. The effective size was 39,159,782 non-N bases. The legacy processing convention documented 200-bp read extension, so this deliberately sparse execution used the labelled `--nomodel --extsize 200` override. Production remains normal single-end MACS3 model building, which the synthetic test executes successfully.

Treatment/input raw reads were 500,000/500,000; fastp retained 391,624/451,279. Mapped reads were 56,459/76,780; MAPQ-qualified evidence was 8,908/10,300; duplicates were 904/2,178; blacklist removal was 209/488; usable evidence was 8,565/9,441. The compact reference lacks chrM, so mitochondrial accounting is zero with `mitochondrial_reference_available = 0`, not an inferred zero fraction from a mitochondrial-containing reference. MACS3 completed with 45 peaks; strict-overlap merging produced 44 consensus regions; FRiP was 0.05149. The signal products and complete count/provenance/status contracts were produced. Differential analysis was disabled.

The run completed in 6 minutes 42 seconds. The largest observed process peak RSS was 1.3 GB (fastp); MACS3 peak RSS was 45.1 MB. Retained H3K27ac fixtures and results occupy about 0.21 GB; total free space after validation remained about 90.4 GB decimal. These values are execution diagnostics only and must not be interpreted as genome-wide mapping quality, ENCODE QC compliance, differential enrichment, enhancer activation, or regulatory mechanism.

## Tool versions

Executed versions were Nextflow 26.04.6, FastQC 0.12.1, fastp 0.23.4, Bowtie2 2.5.5, samtools 1.20, BEDTools 2.31.1, MACS3 3.0.5, deepTools 3.5.6, R 4.5.2, Bioconductor 3.22, DiffBind 3.20.0, DESeq2 1.50.2, MultiQC 1.35, and Python 3.13.13 for deterministic contracts. Container references and immutable digests are recorded in `assets/tool_containers.tsv` and run provenance.
