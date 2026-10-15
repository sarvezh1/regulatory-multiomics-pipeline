process VALIDATE_WGBS_SAMPLESHEET {
    tag 'WGBS samplesheet'
    input: path samplesheet
    output:
    path 'wgbs.samplesheet.validated.tsv', emit: validated
    path 'wgbs.samplesheet.validation.json', emit: report
    path 'raw_inputs.present', optional: true, emit: raw_marker
    path 'processed_inputs.present', optional: true, emit: processed_marker
    path 'validate_wgbs_samplesheet.version.txt', emit: version
    script:
    """
    validate_wgbs_samplesheet.py --input ${samplesheet} --output wgbs.samplesheet.validated.tsv --report wgbs.samplesheet.validation.json
    python --version > validate_wgbs_samplesheet.version.txt 2>&1
    """
}
