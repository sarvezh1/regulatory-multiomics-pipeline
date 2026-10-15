include { VALIDATE_WGBS_SAMPLESHEET } from '../../modules/local/wgbs_validate'
include { VALIDATE_REFERENCE } from '../../modules/local/validate_reference'
include { BISMARK_REFERENCE; WGBS_REFERENCE_CHECKSUM } from '../../modules/local/wgbs_reference'
include { WGBS_FASTQC_RAW; WGBS_FASTQC_TRIMMED } from '../../modules/local/wgbs_fastqc'
include { WGBS_TRIM_GALORE } from '../../modules/local/wgbs_trim'
include { BISMARK_ALIGN } from '../../modules/local/wgbs_align'
include { BISMARK_DEDUP } from '../../modules/local/wgbs_dedup'
include { BISMARK_EXTRACT } from '../../modules/local/wgbs_extract'
include { WGBS_LIBRARY_CONTRACT; ADAPT_PROCESSED_WGBS; AGGREGATE_WGBS_SAMPLE; BUILD_WGBS_UNIVERSE } from '../../modules/local/wgbs_contract'
include { WGBS_QC } from '../../modules/local/wgbs_qc'
include { WGBS_DMRSEQ } from '../../modules/local/wgbs_dm'
include { WGBS_PROCESSING_STATUS } from '../../modules/local/wgbs_status'
include { MULTIQC_WGBS } from '../../modules/local/wgbs_multiqc'
include { BUILD_WGBS_PROVENANCE } from '../../modules/local/wgbs_provenance'

workflow WGBS_WORKFLOW {
    take:
    samplesheet
    fasta
    gtf
    design
    contrast
    run_dm
    route

    main:
    VALIDATE_WGBS_SAMPLESHEET(samplesheet)
    VALIDATE_REFERENCE(fasta, gtf, params.reference_scope)
    WGBS_REFERENCE_CHECKSUM(fasta)
    rows = VALIDATE_WGBS_SAMPLESHEET.out.validated.splitCsv(header:true,sep:'\t')

    if (route == 'raw') {
        raw_reads = rows.filter{it.input_route=='raw'}.map{row->tuple(row as Map,[file(row.read1),file(row.read2)])}
        BISMARK_REFERENCE(fasta, VALIDATE_WGBS_SAMPLESHEET.out.raw_marker)
        WGBS_FASTQC_RAW(raw_reads)
        WGBS_TRIM_GALORE(raw_reads)
        WGBS_FASTQC_TRIMMED(WGBS_TRIM_GALORE.out.trimmed)
        BISMARK_ALIGN(WGBS_TRIM_GALORE.out.trimmed, BISMARK_REFERENCE.out.index)
        BISMARK_DEDUP(BISMARK_ALIGN.out.alignment)
        BISMARK_EXTRACT(BISMARK_DEDUP.out.deduplicated)
        WGBS_LIBRARY_CONTRACT(BISMARK_EXTRACT.out.extraction, fasta)
        WGBS_QC(WGBS_LIBRARY_CONTRACT.out.contract)
        library_contracts = WGBS_LIBRARY_CONTRACT.out.contract.map{meta,contract,rejected,mbias,split,chg,chh,dedup,align,trim->tuple(meta,contract)}
        qc_inputs = WGBS_FASTQC_RAW.out.reports.map{m,z,h->[z,h]}
            .mix(WGBS_FASTQC_TRIMMED.out.reports.map{m,z,h->[z,h]})
            .mix(WGBS_TRIM_GALORE.out.trimmed.map{m,r1,r2,tr->[tr]})
            .mix(BISMARK_ALIGN.out.alignment.map{m,b,a,t->[a]})
            .mix(BISMARK_DEDUP.out.deduplicated.map{m,b,d,a,t->[d]})
            .mix(BISMARK_EXTRACT.out.extraction.map{m,c,mb,s,cg,ch,d,a,t->[mb,s]})
            .flatten().collect()
        route_versions = BISMARK_REFERENCE.out.version.mix(WGBS_FASTQC_RAW.out.version).mix(WGBS_TRIM_GALORE.out.version).mix(WGBS_FASTQC_TRIMMED.out.version).mix(BISMARK_ALIGN.out.version).mix(BISMARK_DEDUP.out.version).mix(BISMARK_EXTRACT.out.version).mix(WGBS_LIBRARY_CONTRACT.out.version).mix(WGBS_QC.out.version)
        reference_checksums = BISMARK_REFERENCE.out.checksums
    } else {
        processed = rows.filter{it.input_route=='processed'}.map{row->tuple(row as Map,file(row.processed_input))}
        ADAPT_PROCESSED_WGBS(processed, fasta)
        library_contracts = ADAPT_PROCESSED_WGBS.out.contract
        qc_inputs = VALIDATE_WGBS_SAMPLESHEET.out.report.collect()
        route_versions = ADAPT_PROCESSED_WGBS.out.version
        reference_checksums = WGBS_REFERENCE_CHECKSUM.out.checksums
    }

    grouped = library_contracts.map{meta,contract->tuple(meta.sample_id,contract)}.groupTuple(by:0)
    AGGREGATE_WGBS_SAMPLE(grouped)
    sample_contracts = AGGREGATE_WGBS_SAMPLE.out.samples.map{sid,contract->contract}.collect()
    BUILD_WGBS_UNIVERSE(sample_contracts, VALIDATE_WGBS_SAMPLESHEET.out.validated)
    WGBS_DMRSEQ(BUILD_WGBS_UNIVERSE.out.long, VALIDATE_WGBS_SAMPLESHEET.out.validated, design, contrast, run_dm)
    WGBS_PROCESSING_STATUS(route, VALIDATE_WGBS_SAMPLESHEET.out.report, sample_contracts)
    MULTIQC_WGBS(qc_inputs, file("${projectDir}/assets/multiqc_config.yml"))
    versions = VALIDATE_WGBS_SAMPLESHEET.out.version.mix(VALIDATE_REFERENCE.out.version).mix(WGBS_REFERENCE_CHECKSUM.out.version).mix(route_versions).mix(AGGREGATE_WGBS_SAMPLE.out.version).mix(BUILD_WGBS_UNIVERSE.out.version).mix(WGBS_DMRSEQ.out.version).mix(WGBS_PROCESSING_STATUS.out.version).mix(MULTIQC_WGBS.out.version).collect()
    BUILD_WGBS_PROVENANCE(VALIDATE_WGBS_SAMPLESHEET.out.validated, VALIDATE_WGBS_SAMPLESHEET.out.report, VALIDATE_REFERENCE.out.report, reference_checksums, WGBS_PROCESSING_STATUS.out.status, WGBS_DMRSEQ.out.status, versions, file("${projectDir}/assets/tool_containers.tsv"), workflow.nextflow.version, workflow.profile ?: 'standard', design, contrast)
}
