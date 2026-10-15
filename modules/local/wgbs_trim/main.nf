process WGBS_TRIM_GALORE {
    tag "${meta.library_id}"
    publishDir "${params.outdir}/wgbs/trimmed", mode: 'copy', overwrite: true, pattern: '*trimming_report.txt'
    input: tuple val(meta), path(reads)
    output:
    tuple val(meta), path('*_val_1.fq.gz'), path('*_val_2.fq.gz'), path('*trimming_report.txt'), emit: trimmed
    path 'trim_galore.version.txt', emit: version
    script:
    """
    trim_galore --paired --quality 20 --length 20 --cores ${task.cpus} --gzip --output_dir . ${reads[0]} ${reads[1]}
    trim_galore --version > trim_galore.version.txt 2>&1
    """
}
