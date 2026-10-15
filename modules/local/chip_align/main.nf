process CHIP_BOWTIE2_ALIGN {
 tag "${meta.library_id}"; publishDir "${params.outdir}/h3k27ac/qc/bowtie2", mode:'copy', overwrite:true, pattern:'*.bowtie2.log'
 input: tuple val(meta), path(reads); path index
 output: tuple val(meta), path("${meta.library_id}.sam"), path("${meta.library_id}.bowtie2.log"), emit:alignment; path 'chip_bowtie2.version.txt', emit:version
 script:
 def input = meta.library_layout == 'paired' ? "--no-mixed --no-discordant -X ${params.chip_paired_max_fragment} -1 ${reads[0]} -2 ${reads[1]}" : "-U ${reads[0]}"
 """
 bowtie2 --very-sensitive --threads ${task.cpus} -x ${index}/reference ${input} -S ${meta.library_id}.sam 2> ${meta.library_id}.bowtie2.log
 bowtie2 --version | head -1 > chip_bowtie2.version.txt
 """
}
