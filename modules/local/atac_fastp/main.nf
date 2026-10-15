process ATAC_FASTP {
    tag "${meta.library_id}"
    publishDir "${params.outdir}/atac/trimmed", mode: 'copy', overwrite: true
    publishDir "${params.outdir}/atac/qc/fastp", mode: 'copy', overwrite: true, pattern: '*.{json,html}'
    input: tuple val(meta), path(reads)
    output: tuple val(meta), path("${meta.library_id}.trimmed_R*.fastq.gz"), path("${meta.library_id}.fastp.json"), path("${meta.library_id}.fastp.html"), emit: trimmed; path 'atac_fastp.version.txt', emit: version
    script:
    """
    fastp --thread ${task.cpus} --in1 ${reads[0]} --in2 ${reads[1]} --out1 ${meta.library_id}.trimmed_R1.fastq.gz --out2 ${meta.library_id}.trimmed_R2.fastq.gz --json ${meta.library_id}.fastp.json --html ${meta.library_id}.fastp.html
    fastp --version > atac_fastp.version.txt 2>&1
    """
}
