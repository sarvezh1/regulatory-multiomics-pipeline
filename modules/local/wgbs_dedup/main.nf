process BISMARK_DEDUP {
    tag "${meta.library_id}"
    publishDir "${params.outdir}/wgbs/alignment/deduplicated", mode: 'copy', overwrite: true
    input: tuple val(meta), path(bam), path(alignment_report), path(trim_report)
    output:
    tuple val(meta), path('*.deduplicated.bam'), path('*deduplication_report.txt'), path(alignment_report), path(trim_report), emit: deduplicated
    path 'bismark_dedup.version.txt', emit: version
    script:
    """
    bismark dedup --paired --bam --parallel ${task.cpus} --output_dir . --outfile ${meta.library_id}.deduplicated ${bam}
    bismark --version > bismark_dedup.version.txt 2>&1
    """
}
