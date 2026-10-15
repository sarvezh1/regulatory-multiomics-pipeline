process BISMARK_EXTRACT {
    tag "${meta.library_id}"
    publishDir "${params.outdir}/wgbs/extraction", mode: 'copy', overwrite: true
    input: tuple val(meta), path(bam), path(dedup_report), path(alignment_report), path(trim_report)
    output:
    tuple val(meta), path('*.bismark.cov.gz'), path('*M-bias.txt'), path('*splitting_report.txt'), path('CHG_context*.txt.gz'), path('CHH_context*.txt.gz'), path(dedup_report), path(alignment_report), path(trim_report), emit: extraction
    path 'bismark_extract.version.txt', emit: version
    script:
    """
    bismark extract --paired-end --gzip --comprehensive --bedGraph --counts --no_overlap --ignore ${params.wgbs_ignore} --ignore_r2 ${params.wgbs_ignore_r2} --ignore_3prime ${params.wgbs_ignore_3prime} --ignore_3prime_r2 ${params.wgbs_ignore_3prime_r2} --output_dir . ${bam}
    bismark --version > bismark_extract.version.txt 2>&1
    """
}
