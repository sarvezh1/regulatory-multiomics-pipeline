process ATAC_PROCESSING_STATUS {
    tag 'ATAC component statuses'
    input: val route; path peaks
    output: path 'atac.processing.status.json', emit: status; path 'atac_status.version.txt', emit: version
    script:
    """
    build_atac_status.py --route '${route}' --peaks ${peaks.join(' ')} --output atac.processing.status.json
    python --version > atac_status.version.txt 2>&1
    """
}
