process ATAC_FASTQC_RAW {
    tag "${meta.library_id}"
    publishDir "${params.outdir}/atac/qc/fastqc/raw", mode: 'copy', overwrite: true
    input: tuple val(meta), path(reads)
    output: tuple val(meta), path('*_fastqc.zip'), path('*_fastqc.html'), emit: reports; path 'atac_fastqc_raw.version.txt', emit: version
    script:
    """
    fastqc --threads ${task.cpus} --outdir . ${reads.join(' ')}
    fastqc --version > atac_fastqc_raw.version.txt 2>&1
    """
}
process ATAC_FASTQC_TRIMMED {
    tag "${meta.library_id}"
    publishDir "${params.outdir}/atac/qc/fastqc/trimmed", mode: 'copy', overwrite: true
    input: tuple val(meta), path(reads), path(json), path(html)
    output: tuple val(meta), path('*_fastqc.zip'), path('*_fastqc.html'), emit: reports; path 'atac_fastqc_trimmed.version.txt', emit: version
    script:
    """
    fastqc --threads ${task.cpus} --outdir . ${reads.join(' ')}
    fastqc --version > atac_fastqc_trimmed.version.txt 2>&1
    """
}
