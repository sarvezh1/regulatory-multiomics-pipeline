process ATAC_MERGE_SAMPLE {
    tag "${sample_id}"
    publishDir "${params.outdir}/atac/alignment/samples", mode: 'copy', overwrite: true, pattern: '*.{bam,bai}'
    input: tuple val(sample_id), val(metas), path(bams), path(bais), path(fragments), path(blacklist_metrics), path(filter_metrics)
    output: tuple val(sample_id), val(metas), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), path("${sample_id}.fragments.bed"), path(blacklist_metrics), path(filter_metrics), emit: samples; path 'samtools_merge.version.txt', emit: version
    script:
    def merge = bams.size() == 1 ? "cp ${bams[0]} ${sample_id}.bam" : "samtools merge -@ ${task.cpus} -f ${sample_id}.bam ${bams.join(' ')}"
    """
    ${merge}
    samtools index ${sample_id}.bam
    cat ${fragments.join(' ')} | sort -k1,1 -k2,2n -k3,3n -k4,4 > ${sample_id}.fragments.bed
    samtools --version | head -1 > samtools_merge.version.txt
    """
}
