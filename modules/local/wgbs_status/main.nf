process WGBS_PROCESSING_STATUS {
    tag 'WGBS component statuses'
    publishDir "${params.outdir}/wgbs/provenance", mode: 'copy', overwrite: true
    input: val route; path validation_report; path contracts
    output: path 'wgbs.processing_status.json', emit: status; path 'wgbs_status.version.txt', emit: version
    script:
    """
    build_wgbs_status.py --route ${route} --validation-report ${validation_report} --contracts ${contracts.join(' ')} --output wgbs.processing_status.json
    python --version > wgbs_status.version.txt 2>&1
    """
}
