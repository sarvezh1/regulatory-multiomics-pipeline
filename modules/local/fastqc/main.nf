process FASTQC_RAW {
    tag "${meta.library_id}"
    publishDir "${params.outdir}/rna/qc/fastqc/raw", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*_fastqc.zip'), path('*_fastqc.html'), emit: reports
    path 'fastqc.version.txt', emit: version

    script:
    """
    fastqc --threads ${task.cpus} --outdir . ${reads.join(' ')}
    fastqc --version > fastqc.version.txt 2>&1
    """
}

process FASTQC_TRIMMED {
    tag "${meta.library_id}"
    publishDir "${params.outdir}/rna/qc/fastqc/trimmed", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(reads), path(fastp_json), path(fastp_html)

    output:
    tuple val(meta), path('*_fastqc.zip'), path('*_fastqc.html'), emit: reports
    path 'fastqc.version.txt', emit: version

    script:
    """
    fastqc --threads ${task.cpus} --outdir . ${reads.join(' ')}
    fastqc --version > fastqc.version.txt 2>&1
    """
}
