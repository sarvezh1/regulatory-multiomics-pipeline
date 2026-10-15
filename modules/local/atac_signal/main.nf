process ATAC_BIGWIG {
    tag "${sample_id}"
    publishDir "${params.outdir}/atac/signal", mode: 'copy', overwrite: true
    input: tuple val(sample_id), val(metas), path(bam), path(bai), path(fragments), path(peaks), path(summits), path(xls), path(log), path(blacklist_metrics), path(filter_metrics)
    output: tuple val(sample_id), path("${sample_id}.fragment_coverage.CPM.bw"), emit: signal; path 'deeptools.version.txt', emit: version
    script:
    """
    bamCoverage --bam ${bam} --outFileName ${sample_id}.fragment_coverage.CPM.bw --outFileFormat bigwig --normalizeUsing CPM --extendReads --binSize 10 --numberOfProcessors ${task.cpus}
    bamCoverage --version > deeptools.version.txt
    """
}
