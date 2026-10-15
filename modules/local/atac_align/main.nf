process BOWTIE2_ALIGN {
    tag "${meta.library_id}"
    publishDir "${params.outdir}/atac/qc/bowtie2", mode: 'copy', overwrite: true, pattern: '*.bowtie2.log'
    input: tuple val(meta), path(reads); path index
    output: tuple val(meta), path("${meta.library_id}.sam"), path("${meta.library_id}.bowtie2.log"), emit: alignment; path 'bowtie2.version.txt', emit: version
    script:
    """
    bowtie2 --very-sensitive -X 2000 --no-mixed --no-discordant --threads ${task.cpus} -x ${index}/reference -1 ${reads[0]} -2 ${reads[1]} -S ${meta.library_id}.sam 2> ${meta.library_id}.bowtie2.log
    bowtie2 --version | head -1 > bowtie2.version.txt
    """
}
