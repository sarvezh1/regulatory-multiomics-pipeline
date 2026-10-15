include { VALIDATE_ATAC_SAMPLESHEET } from '../../modules/local/atac_validate'
include { VALIDATE_REFERENCE } from '../../modules/local/validate_reference'
include { BOWTIE2_INDEX; DERIVE_ATAC_TSS } from '../../modules/local/atac_reference'
include { ATAC_FASTQC_RAW; ATAC_FASTQC_TRIMMED } from '../../modules/local/atac_fastqc'
include { ATAC_FASTP } from '../../modules/local/atac_fastp'
include { BOWTIE2_ALIGN } from '../../modules/local/atac_align'
include { ATAC_SAMTOOLS_PREP; ATAC_BLACKLIST_FILTER; ATAC_FINAL_BAM } from '../../modules/local/atac_filter'
include { ATAC_MERGE_SAMPLE } from '../../modules/local/atac_sample'
include { MACS3_ATAC } from '../../modules/local/atac_peaks'
include { ATAC_BIGWIG } from '../../modules/local/atac_signal'
include { ATAQV; ATAC_CUSTOM_QC } from '../../modules/local/atac_qc'
include { BUILD_ATAC_CONSENSUS; BUILD_ATAC_MATRIX; ADAPT_PROCESSED_ATAC } from '../../modules/local/atac_contract'
include { ATAC_DIFFBIND } from '../../modules/local/atac_da'
include { ATAC_PROCESSING_STATUS } from '../../modules/local/atac_status'
include { ATAC_COUNT_REPLICATE_QC; ATAC_PEAK_JACCARD } from '../../modules/local/atac_replicate_qc'
include { MULTIQC_ATAC } from '../../modules/local/atac_multiqc'
include { BUILD_ATAC_PROVENANCE } from '../../modules/local/atac_provenance'

workflow ATAC_WORKFLOW {
    take:
    samplesheet
    fasta
    gtf
    blacklist
    design
    contrast
    run_da
    route

    main:
    VALIDATE_ATAC_SAMPLESHEET(samplesheet)
    VALIDATE_REFERENCE(fasta, gtf, params.reference_scope)
    rows = VALIDATE_ATAC_SAMPLESHEET.out.validated.splitCsv(header:true,sep:'\t')

    if (route == 'raw') {
        raw_reads = rows.filter { it.input_route == 'raw' }.map { row -> tuple(row as Map, [file(row.read1),file(row.read2)]) }
        BOWTIE2_INDEX(fasta, VALIDATE_ATAC_SAMPLESHEET.out.raw_marker)
        DERIVE_ATAC_TSS(gtf, VALIDATE_ATAC_SAMPLESHEET.out.raw_marker)
        ATAC_FASTQC_RAW(raw_reads)
        ATAC_FASTP(raw_reads)
        ATAC_FASTQC_TRIMMED(ATAC_FASTP.out.trimmed)
        trimmed_for_alignment = ATAC_FASTP.out.trimmed.map { meta,reads,json,html -> tuple(meta,reads) }
        BOWTIE2_ALIGN(trimmed_for_alignment, BOWTIE2_INDEX.out.index)
        ATAC_SAMTOOLS_PREP(BOWTIE2_ALIGN.out.alignment)
        ATAC_BLACKLIST_FILTER(ATAC_SAMTOOLS_PREP.out.prepared, blacklist)
        ATAC_FINAL_BAM(ATAC_BLACKLIST_FILTER.out.filtered_names)
        grouped = ATAC_FINAL_BAM.out.filtered.map { meta,bam,bai,frag,bm,fm -> tuple(meta.sample_id,meta,bam,bai,frag,bm,fm) }.groupTuple(by:0)
        ATAC_MERGE_SAMPLE(grouped)
        MACS3_ATAC(ATAC_MERGE_SAMPLE.out.samples)
        ATAC_BIGWIG(MACS3_ATAC.out.peaks)
        ATAQV(MACS3_ATAC.out.peaks, DERIVE_ATAC_TSS.out.bed)
        ATAC_CUSTOM_QC(MACS3_ATAC.out.peaks, DERIVE_ATAC_TSS.out.bed)
        peak_files = MACS3_ATAC.out.peaks.map { sid,metas,bam,bai,frag,peaks,summits,xls,log,bm,fm -> peaks }.collect()
        fragment_files = MACS3_ATAC.out.peaks.map { sid,metas,bam,bai,frag,peaks,summits,xls,log,bm,fm -> frag }.collect()
        BUILD_ATAC_CONSENSUS(peak_files)
        BUILD_ATAC_MATRIX(BUILD_ATAC_CONSENSUS.out.regions, fragment_files)
        counts_ch = BUILD_ATAC_MATRIX.out.matrix
        regions_ch = BUILD_ATAC_CONSENSUS.out.regions
        ATAC_PROCESSING_STATUS(route, peak_files)
        ATAC_PEAK_JACCARD(peak_files)
        tss_for_provenance = DERIVE_ATAC_TSS.out.bed
        qc_inputs = ATAC_FASTQC_RAW.out.reports.map{m,z,h->[z,h]}
            .mix(ATAC_FASTQC_TRIMMED.out.reports.map{m,z,h->[z,h]})
            .mix(ATAC_FASTP.out.trimmed.map{m,r,j,h->[j,h]})
            .mix(BOWTIE2_ALIGN.out.alignment.map{m,s,l->[l]})
            .mix(ATAC_SAMTOOLS_PREP.out.qc.map{m,f,s->[f,s]})
            .flatten().collect()
        raw_versions = BOWTIE2_INDEX.out.version.mix(DERIVE_ATAC_TSS.out.version).mix(ATAC_FASTQC_RAW.out.version).mix(ATAC_FASTP.out.version).mix(ATAC_FASTQC_TRIMMED.out.version).mix(BOWTIE2_ALIGN.out.version).mix(ATAC_SAMTOOLS_PREP.out.version).mix(ATAC_BLACKLIST_FILTER.out.version).mix(ATAC_FINAL_BAM.out.version).mix(ATAC_MERGE_SAMPLE.out.version).mix(MACS3_ATAC.out.version).mix(ATAC_BIGWIG.out.version).mix(ATAQV.out.version).mix(ATAC_CUSTOM_QC.out.version).mix(BUILD_ATAC_CONSENSUS.out.version).mix(BUILD_ATAC_MATRIX.out.version).mix(ATAC_PEAK_JACCARD.out.version)
    } else {
        proc_matrix = rows.filter{it.input_route=='processed'}.map{row->file(row.processed_input)}.first()
        ADAPT_PROCESSED_ATAC(proc_matrix, VALIDATE_ATAC_SAMPLESHEET.out.validated)
        counts_ch = ADAPT_PROCESSED_ATAC.out.counts
        regions_ch = ADAPT_PROCESSED_ATAC.out.regions
        ATAC_PROCESSING_STATUS(route, blacklist)
        tss_for_provenance = gtf
        qc_inputs = VALIDATE_ATAC_SAMPLESHEET.out.report.collect()
        raw_versions = ADAPT_PROCESSED_ATAC.out.version
    }

    sample_metadata = VALIDATE_ATAC_SAMPLESHEET.out.validated
    ATAC_COUNT_REPLICATE_QC(counts_ch)
    ATAC_DIFFBIND(counts_ch, sample_metadata, design, contrast, run_da)
    MULTIQC_ATAC(qc_inputs, file("${projectDir}/assets/multiqc_config.yml"))
    versions = VALIDATE_ATAC_SAMPLESHEET.out.version.mix(VALIDATE_REFERENCE.out.version).mix(raw_versions).mix(ATAC_COUNT_REPLICATE_QC.out.version).mix(ATAC_DIFFBIND.out.version).mix(ATAC_PROCESSING_STATUS.out.version).mix(MULTIQC_ATAC.out.version).collect()
    BUILD_ATAC_PROVENANCE(VALIDATE_ATAC_SAMPLESHEET.out.validated, VALIDATE_ATAC_SAMPLESHEET.out.report, VALIDATE_REFERENCE.out.report, blacklist, tss_for_provenance, ATAC_PROCESSING_STATUS.out.status, ATAC_DIFFBIND.out.status, versions, file("${projectDir}/assets/tool_containers.tsv"), workflow.nextflow.version, workflow.profile ?: 'standard', design, contrast, params.atac_effective_genome_size)
}
