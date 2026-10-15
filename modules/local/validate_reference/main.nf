process VALIDATE_REFERENCE {
    tag 'reference'

    input:
    path fasta
    path gtf
    val reference_scope

    output:
    path 'reference.validation.json', emit: report
    path 'reference.validated', emit: marker
    path 'validate_reference.version.txt', emit: version

    script:
    """
    validate_reference.py --fasta ${fasta} --gtf ${gtf} --scope '${reference_scope}' --output reference.validation.json
    touch reference.validated
    python --version > validate_reference.version.txt 2>&1
    """
}
