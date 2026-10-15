process SAMTOOLS_QC {
    tag "${meta.library_id}"
    publishDir "${params.outdir}/rna/alignment", mode: 'copy', overwrite: true, pattern: '*.{bai,flagstat.txt,stats.txt}'

    input:
    tuple val(meta), path(bam), path(log_final), path(log_out), path(sj)

    output:
    tuple val(meta), path("${bam}.bai"), path("${meta.library_id}.flagstat.txt"), path("${meta.library_id}.stats.txt"), emit: qc
    path 'samtools.version.txt', emit: version

    script:
    """
    samtools index -@ ${task.cpus} ${bam}
    samtools flagstat -@ ${task.cpus} ${bam} > ${meta.library_id}.flagstat.txt
    samtools stats -@ ${task.cpus} ${bam} > ${meta.library_id}.stats.txt
    samtools --version | head -n 1 > samtools.version.txt
    """
}
