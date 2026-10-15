include { VALIDATE_REFERENCE } from '../../modules/local/validate_reference'
include { BUILD_INTEGRATION_REFERENCE } from '../../modules/local/integration_reference'
include { RUN_REGULATORY_INTEGRATION } from '../../modules/local/integration_run'

workflow INTEGRATION_WORKFLOW {
    take:
    manifest
    fasta
    gtf
    mode
    match_type
    fdr
    distal_max

    main:
    VALIDATE_REFERENCE(fasta, gtf, params.reference_scope)
    BUILD_INTEGRATION_REFERENCE(gtf, fasta, params.annotation_name)
    contract_inputs = Channel.of(manifest).splitCsv(header:true, sep:'\t').map { row -> file(row.path, checkIfExists: true) }.collect()
    RUN_REGULATORY_INTEGRATION(manifest, contract_inputs, BUILD_INTEGRATION_REFERENCE.out.reference, mode, match_type, fdr, distal_max)

    emit:
    results = RUN_REGULATORY_INTEGRATION.out.results
}
