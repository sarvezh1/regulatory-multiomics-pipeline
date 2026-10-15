process CHIP_FASTQC_RAW {
 tag "${meta.library_id}"; publishDir "${params.outdir}/h3k27ac/qc/fastqc/raw", mode:'copy', overwrite:true
 input: tuple val(meta), path(reads)
 output: tuple val(meta), path('*_fastqc.zip'), path('*_fastqc.html'), emit:reports; path 'chip_fastqc_raw.version.txt', emit:version
 script:
 """
 fastqc --threads ${task.cpus} --outdir . ${reads.join(' ')}
 fastqc --version > chip_fastqc_raw.version.txt 2>&1
 """
}
process CHIP_FASTQC_TRIMMED {
 tag "${meta.library_id}"; publishDir "${params.outdir}/h3k27ac/qc/fastqc/trimmed", mode:'copy', overwrite:true
 input: tuple val(meta), path(reads), path(json), path(html)
 output: tuple val(meta), path('*_fastqc.zip'), path('*_fastqc.html'), emit:reports; path 'chip_fastqc_trimmed.version.txt', emit:version
 script:
 """
 fastqc --threads ${task.cpus} --outdir . ${reads.join(' ')}
 fastqc --version > chip_fastqc_trimmed.version.txt 2>&1
 """
}
