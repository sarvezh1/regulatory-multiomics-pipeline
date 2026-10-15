process VALIDATE_ATAC_SAMPLESHEET {
    tag 'ATAC samplesheet'
    input: path samplesheet
    output:
    path 'atac.samplesheet.validated.tsv', emit: validated
    path 'atac.samplesheet.validation.json', emit: report
    path 'raw_inputs.present', optional: true, emit: raw_marker
    path 'processed_inputs.present', optional: true, emit: processed_marker
    path 'validate_atac_samplesheet.version.txt', emit: version
    script:
    """
    validate_atac_samplesheet.py --input ${samplesheet} --output atac.samplesheet.validated.tsv --report atac.samplesheet.validation.json
    python --version > validate_atac_samplesheet.version.txt 2>&1
    """
}
