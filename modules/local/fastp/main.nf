process FASTP {
    tag "${meta.library_id}"
    publishDir "${params.outdir}/rna/trimmed", mode: 'copy', overwrite: true
    publishDir "${params.outdir}/rna/qc/fastp", mode: 'copy', overwrite: true, pattern: '*.{json,html}'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.library_id}.trimmed*.fastq.gz"), path("${meta.library_id}.fastp.json"), path("${meta.library_id}.fastp.html"), emit: trimmed
    path 'fastp.version.txt', emit: version

    script:
    def command = meta.library_layout == 'paired' ?
        "fastp --thread ${task.cpus} --in1 ${reads[0]} --in2 ${reads[1]} --out1 ${meta.library_id}.trimmed_R1.fastq.gz --out2 ${meta.library_id}.trimmed_R2.fastq.gz --json ${meta.library_id}.fastp.json --html ${meta.library_id}.fastp.html" :
        "fastp --thread ${task.cpus} --in1 ${reads[0]} --out1 ${meta.library_id}.trimmed.fastq.gz --json ${meta.library_id}.fastp.json --html ${meta.library_id}.fastp.html"
    """
    ${command}
    fastp --version > fastp.version.txt 2>&1
    """
}
