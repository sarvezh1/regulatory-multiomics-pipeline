process WGBS_FASTQC_RAW {
    tag "${meta.library_id}"
    publishDir "${params.outdir}/wgbs/qc/fastqc/raw", mode: 'copy', overwrite: true
    input: tuple val(meta), path(reads)
    output: tuple val(meta), path('*_fastqc.zip'), path('*_fastqc.html'), emit: reports
    path 'wgbs_fastqc_raw.version.txt', emit: version
    script:
    """
    fastqc --threads ${task.cpus} --quiet ${reads.join(' ')}
    fastqc --version > wgbs_fastqc_raw.version.txt 2>&1
    """
}

process WGBS_FASTQC_TRIMMED {
    tag "${meta.library_id}"
    publishDir "${params.outdir}/wgbs/qc/fastqc/trimmed", mode: 'copy', overwrite: true
    input: tuple val(meta), path(read1), path(read2), path(trim_report)
    output: tuple val(meta), path('*_fastqc.zip'), path('*_fastqc.html'), emit: reports
    path 'wgbs_fastqc_trimmed.version.txt', emit: version
    script:
    """
    fastqc --threads ${task.cpus} --quiet ${read1} ${read2}
    fastqc --version > wgbs_fastqc_trimmed.version.txt 2>&1
    """
}
