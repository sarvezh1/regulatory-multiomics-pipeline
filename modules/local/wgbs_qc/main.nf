process WGBS_QC {
    tag "${meta.library_id}"
    publishDir "${params.outdir}/wgbs/qc/metrics", mode: 'copy', overwrite: true
    input: tuple val(meta), path(contract), path(rejected), path(mbias), path(split_report), path(chg), path(chh), path(dedup_report), path(alignment_report), path(trim_report)
    output:
    tuple val(meta), path("${meta.library_id}.wgbs_qc.tsv"), path("${meta.library_id}.wgbs_qc.json"), emit: metrics
    path 'wgbs_qc.version.txt', emit: version
    script:
    def control = meta.conversion_control_contig ?: ''
    def trimArgs = trim_report.collect { "--trim-report ${it}" }.join(' ')
    """
    wgbs_qc.py --contract ${contract} --rejected ${rejected} --mbias ${mbias} --split-report ${split_report} --chg ${chg} --chh ${chh} --dedup-report ${dedup_report} --alignment-report ${alignment_report} ${trimArgs} --sample-id ${meta.sample_id} --library-id ${meta.library_id} --min-coverage ${params.wgbs_min_coverage} --conversion-control-contig '${control}' --output ${meta.library_id}.wgbs_qc.tsv --json-output ${meta.library_id}.wgbs_qc.json
    python --version > wgbs_qc.version.txt 2>&1
    """
}
