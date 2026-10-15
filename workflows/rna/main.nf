include { VALIDATE_SAMPLESHEET } from '../../modules/local/validate_samplesheet'
include { VALIDATE_REFERENCE } from '../../modules/local/validate_reference'
include { STAR_INDEX } from '../../modules/local/star_index'
include { FASTQC_RAW; FASTQC_TRIMMED } from '../../modules/local/fastqc'
include { FASTP } from '../../modules/local/fastp'
include { STAR_ALIGN } from '../../modules/local/star_align'
include { SAMTOOLS_QC } from '../../modules/local/samtools_qc'
include { FEATURECOUNTS } from '../../modules/local/featurecounts'
include { FEATURECOUNTS_TO_CONTRACT; ADAPT_PROCESSED_COUNTS; MERGE_RNA_COUNTS } from '../../modules/local/rna_contract'
include { RNA_DESEQ2 } from '../../modules/local/deseq2'
include { MULTIQC_RNA } from '../../modules/local/multiqc'
include { BUILD_PROVENANCE } from '../../modules/local/provenance'

workflow RNA_WORKFLOW {
    take:
    samplesheet
    fasta
    gtf
    design_formula
    contrast
    run_de

    main:
    VALIDATE_SAMPLESHEET(samplesheet)
    VALIDATE_REFERENCE(fasta, gtf, params.reference_scope)

    rows = VALIDATE_SAMPLESHEET.out.validated
        .splitCsv(header: true, sep: '\t')

    raw_reads = rows
        .filter { it.input_route == 'raw' }
        .map { row ->
            def meta = row as Map
            def reads = row.library_layout == 'paired' ? [file(row.read1), file(row.read2)] : [file(row.read1)]
            tuple(meta, reads)
        }

    processed_counts = rows
        .filter { it.input_route == 'processed' }
        .map { row -> tuple(row as Map, file(row.processed_input)) }

    STAR_INDEX(fasta, gtf, VALIDATE_SAMPLESHEET.out.raw_marker)
    FASTQC_RAW(raw_reads)
    FASTP(raw_reads)
    FASTQC_TRIMMED(FASTP.out.trimmed)
    STAR_ALIGN(FASTP.out.trimmed, STAR_INDEX.out.index)
    SAMTOOLS_QC(STAR_ALIGN.out.alignment)
    FEATURECOUNTS(STAR_ALIGN.out.alignment, gtf)
    FEATURECOUNTS_TO_CONTRACT(FEATURECOUNTS.out.counts, gtf)
    ADAPT_PROCESSED_COUNTS(processed_counts)

    contracts = FEATURECOUNTS_TO_CONTRACT.out.contract
        .mix(ADAPT_PROCESSED_COUNTS.out.contract)
        .map { meta, contract -> contract }
        .collect()
    MERGE_RNA_COUNTS(contracts, VALIDATE_SAMPLESHEET.out.validated)

    RNA_DESEQ2(
        MERGE_RNA_COUNTS.out.matrix,
        MERGE_RNA_COUNTS.out.metadata,
        MERGE_RNA_COUNTS.out.annotations,
        design_formula,
        contrast,
        run_de
    )

    qc_inputs = FASTQC_RAW.out.reports
        .map { meta, zip, html -> [zip, html] }
        .mix(FASTQC_TRIMMED.out.reports.map { meta, zip, html -> [zip, html] })
        .mix(FASTP.out.trimmed.map { meta, reads, json, html -> [json, html] })
        .mix(STAR_ALIGN.out.alignment.map { meta, bam, finalLog, log, sj -> [finalLog] })
        .mix(SAMTOOLS_QC.out.qc.map { meta, bai, flagstat, stats -> [flagstat, stats] })
        .mix(FEATURECOUNTS.out.counts.map { meta, counts, summary -> [summary] })
        .mix(Channel.value(file("${projectDir}/assets/rna_contract_mqc.yaml")))
        .flatten()
        .collect()
    MULTIQC_RNA(qc_inputs, file("${projectDir}/assets/multiqc_config.yml"))

    version_files = VALIDATE_SAMPLESHEET.out.version
        .mix(VALIDATE_REFERENCE.out.version)
        .mix(STAR_INDEX.out.version)
        .mix(FASTQC_RAW.out.version)
        .mix(FASTP.out.version)
        .mix(FASTQC_TRIMMED.out.version)
        .mix(STAR_ALIGN.out.version)
        .mix(SAMTOOLS_QC.out.version)
        .mix(FEATURECOUNTS.out.version)
        .mix(FEATURECOUNTS_TO_CONTRACT.out.version)
        .mix(ADAPT_PROCESSED_COUNTS.out.version)
        .mix(MERGE_RNA_COUNTS.out.version)
        .mix(RNA_DESEQ2.out.version)
        .mix(MULTIQC_RNA.out.version)
        .collect()

    BUILD_PROVENANCE(
        VALIDATE_SAMPLESHEET.out.validated,
        VALIDATE_REFERENCE.out.report,
        RNA_DESEQ2.out.status,
        RNA_DESEQ2.out.design,
        version_files,
        file("${projectDir}/assets/tool_containers.tsv"),
        workflow.nextflow.version,
        design_formula,
        contrast,
        run_de,
        workflow.profile ?: 'standard'
    )
}
