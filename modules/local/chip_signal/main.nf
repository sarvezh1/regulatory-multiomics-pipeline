process CHIP_SIGNAL {
 tag "${treatment_id}"; publishDir "${params.outdir}/h3k27ac/signal", mode:'copy', overwrite:true
 input: tuple val(treatment_id), val(tm), path(tb), path(tbi), path(te), path(tmet), val(control_id), val(cm), path(cb), path(cbi), path(ce), path(cmet), path(peaks), path(summits), path(xls), path(pileup), path(lambda), path(fe), path(log)
 output: tuple path("${treatment_id}.ChIP.CPM.bw"), path("${control_id}.input.CPM.bw"), path("${treatment_id}.input_normalized.log2.bw"), path(fe), emit:signal; path 'chip_deeptools.version.txt', emit:version
 script:
 """
 bamCoverage -b ${tb} -o ${treatment_id}.ChIP.CPM.bw --normalizeUsing CPM --binSize 25 --numberOfProcessors ${task.cpus}
 bamCoverage -b ${cb} -o ${control_id}.input.CPM.bw --normalizeUsing CPM --binSize 25 --numberOfProcessors ${task.cpus}
 bamCompare -b1 ${tb} -b2 ${cb} -o ${treatment_id}.input_normalized.log2.bw --operation log2 --normalizeUsing CPM --scaleFactorsMethod None --binSize 25 --numberOfProcessors ${task.cpus}
 bamCoverage --version > chip_deeptools.version.txt 2>&1
 """
}
