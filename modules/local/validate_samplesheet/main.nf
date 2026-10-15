process VALIDATE_SAMPLESHEET {
    tag 'samplesheet'

    input:
    path samplesheet

    output:
    path 'samplesheet.validated.tsv', emit: validated
    path 'samplesheet.validation.json', emit: report
    path 'raw_inputs.present', optional: true, emit: raw_marker
    path 'validate_samplesheet.version.txt', emit: version

    script:
    """
    validate_samplesheet.py --input ${samplesheet} --output samplesheet.validated.tsv --report samplesheet.validation.json
    python --version > validate_samplesheet.version.txt 2>&1
    """
}
