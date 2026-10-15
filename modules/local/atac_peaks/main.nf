process MACS3_ATAC {
    tag "${sample_id}"
    publishDir "${params.outdir}/atac/peaks/samples", mode: 'copy', overwrite: true
    input: tuple val(sample_id), val(metas), path(bam), path(bai), path(fragments), path(blacklist_metrics), path(filter_metrics)
    output: tuple val(sample_id), val(metas), path(bam), path(bai), path(fragments), path("${sample_id}_peaks.narrowPeak"), path("${sample_id}_summits.bed"), path("${sample_id}_peaks.xls"), path("${sample_id}.macs3.log"), path(blacklist_metrics), path(filter_metrics), emit: peaks; path 'macs3.version.txt', emit: version
    script:
    """
    macs3 callpeak -t ${bam} -f BAMPE -g ${params.atac_effective_genome_size} -n ${sample_id} -q 0.01 --keep-dup all --call-summits 2> ${sample_id}.macs3.log
    test -f ${sample_id}_peaks.narrowPeak || touch ${sample_id}_peaks.narrowPeak
    test -f ${sample_id}_summits.bed || touch ${sample_id}_summits.bed
    test -f ${sample_id}_peaks.xls || touch ${sample_id}_peaks.xls
    macs3 --version > macs3.version.txt
    """
}
