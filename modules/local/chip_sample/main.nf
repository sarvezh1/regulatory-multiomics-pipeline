process CHIP_MERGE_SAMPLE {
 tag "${sample_id}"; publishDir "${params.outdir}/h3k27ac/alignment/samples", mode:'copy', overwrite:true, pattern:'*.{bam,bai}'
 input: tuple val(sample_id), val(metas), path(bams), path(bais), path(evidence), path(metrics)
 output: tuple val(sample_id), val(metas), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), path("${sample_id}.evidence.bed"), path(metrics), emit:samples; path 'chip_merge.version.txt', emit:version
 script:
 def merge=bams.size()==1 ? "cp ${bams[0]} ${sample_id}.bam" : "samtools merge -f ${sample_id}.bam ${bams.join(' ')}"
 """
 ${merge}
 samtools index ${sample_id}.bam
 cat ${evidence.join(' ')} | sort -k1,1 -k2,2n > ${sample_id}.evidence.bed
 samtools --version | head -1 > chip_merge.version.txt
 """
}
