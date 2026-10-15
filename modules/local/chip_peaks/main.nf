process MACS3_H3K27AC {
 tag "${treatment_id}"; publishDir "${params.outdir}/h3k27ac/peaks/samples", mode:'copy', overwrite:true
 input: tuple val(treatment_id), val(treatment_meta), path(treatment_bam), path(treatment_bai), path(treatment_evidence), path(treatment_metrics), val(control_id), val(control_meta), path(control_bam), path(control_bai), path(control_evidence), path(control_metrics)
 output: tuple val(treatment_id), val(treatment_meta), path(treatment_bam), path(treatment_bai), path(treatment_evidence), path(treatment_metrics), val(control_id), val(control_meta), path(control_bam), path(control_bai), path(control_evidence), path(control_metrics), path("${treatment_id}_peaks.narrowPeak"), path("${treatment_id}_summits.bed"), path("${treatment_id}_peaks.xls"), path("${treatment_id}_treat_pileup.bdg"), path("${treatment_id}_control_lambda.bdg"), path("${treatment_id}.fold_enrichment.bdg"), path("${treatment_id}.macs3.log"), emit:peaks; path 'chip_macs3.version.txt', emit:version
 script:
 def format=treatment_meta[0].library_layout=='paired' ? 'BAMPE' : 'BAM'
 def model=treatment_meta[0].library_layout=='single' && params.chip_se_validation_extsize.toInteger()>0 ? "--nomodel --extsize ${params.chip_se_validation_extsize}" : ''
 """
 macs3 callpeak -t ${treatment_bam} -c ${control_bam} -f ${format} -g ${params.chip_effective_genome_size} -n ${treatment_id} -q 0.01 --keep-dup all --call-summits -B --SPMR ${model} 2> ${treatment_id}.macs3.log
 test -f ${treatment_id}_peaks.narrowPeak || touch ${treatment_id}_peaks.narrowPeak
 test -f ${treatment_id}_summits.bed || touch ${treatment_id}_summits.bed
 test -f ${treatment_id}_peaks.xls || touch ${treatment_id}_peaks.xls
 test -f ${treatment_id}_treat_pileup.bdg || touch ${treatment_id}_treat_pileup.bdg
 test -f ${treatment_id}_control_lambda.bdg || touch ${treatment_id}_control_lambda.bdg
 if [ -s ${treatment_id}_treat_pileup.bdg ] && [ -s ${treatment_id}_control_lambda.bdg ]; then macs3 bdgcmp -t ${treatment_id}_treat_pileup.bdg -c ${treatment_id}_control_lambda.bdg -m FE -o ${treatment_id}.fold_enrichment.bdg; else touch ${treatment_id}.fold_enrichment.bdg; fi
 macs3 --version > chip_macs3.version.txt
 """
}
