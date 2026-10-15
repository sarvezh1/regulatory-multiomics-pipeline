process CHIP_QC {
 tag "${treatment_id}"; publishDir "${params.outdir}/h3k27ac/qc/metrics", mode:'copy', overwrite:true
 input: tuple val(treatment_id), val(tm), path(tb), path(tbi), path(te), path(tmet), val(control_id), val(cm), path(cb), path(cbi), path(ce), path(cmet), path(peaks), path(summits), path(xls), path(pileup), path(lambda), path(fe), path(log)
 output: path("${treatment_id}.chip_qc.json"), emit:metrics; path 'chip_qc.version.txt', emit:version
 script:
 """
 chip_qc.py --sample ${treatment_id} --metrics ${tmet.join(' ')} --peaks ${peaks} --evidence ${te} --output ${treatment_id}.chip_qc.json
 python --version > chip_qc.version.txt 2>&1
 """
}
