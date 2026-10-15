process CHIP_FASTP {
 tag "${meta.library_id}"; publishDir "${params.outdir}/h3k27ac/trimmed", mode:'copy', overwrite:true; publishDir "${params.outdir}/h3k27ac/qc/fastp", mode:'copy', overwrite:true, pattern:'*.{json,html}'
 input: tuple val(meta), path(reads)
 output: tuple val(meta), path("${meta.library_id}.trimmed*.fastq.gz"), path("${meta.library_id}.fastp.json"), path("${meta.library_id}.fastp.html"), emit:trimmed; path 'chip_fastp.version.txt', emit:version
 script:
 def cmd = meta.library_layout == 'paired' ? "fastp --thread ${task.cpus} --in1 ${reads[0]} --in2 ${reads[1]} --out1 ${meta.library_id}.trimmed_R1.fastq.gz --out2 ${meta.library_id}.trimmed_R2.fastq.gz" : "fastp --thread ${task.cpus} --in1 ${reads[0]} --out1 ${meta.library_id}.trimmed.fastq.gz"
 """
 ${cmd} --json ${meta.library_id}.fastp.json --html ${meta.library_id}.fastp.html
 fastp --version > chip_fastp.version.txt 2>&1
 """
}
