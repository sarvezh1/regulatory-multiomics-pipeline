process BUILD_INTEGRATION_REFERENCE {
    tag 'GENCODE v50 TSS and promoters'

    input:
    path gtf
    path fasta
    val annotation

    output:
    path 'integration_reference', emit: reference

    script:
    """
    mkdir -p integration_reference
    build_integration_reference.R ${gtf} ${fasta} integration_reference '${annotation}'
    """
}
