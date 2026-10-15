process CHIP_PROCESSING_STATUS {
 tag 'H3K27ac statuses'
 input: val route; path peaks
 output: path 'h3k27ac.processing.status.json', emit:status; path 'chip_status.version.txt', emit:version
 script:
 """
 build_chip_status.py --route '${route}' --peaks ${peaks.join(' ')} --output h3k27ac.processing.status.json
 python --version > chip_status.version.txt 2>&1
 """
}
