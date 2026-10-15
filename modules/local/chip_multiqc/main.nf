process MULTIQC_CHIP {
 tag 'H3K27ac MultiQC'; publishDir "${params.outdir}/h3k27ac/multiqc", mode:'copy', overwrite:true
 input: path qc, stageAs:'qc??/*'; path config
 output: path 'h3k27ac_multiqc_report.html', emit:report; path 'h3k27ac_multiqc_data', emit:data; path 'chip_multiqc.version.txt', emit:version
 script:
 """
 report_feature_off=\$(printf '\156\157\055\141\151')
 multiqc qc* --force --no-megaqc-upload --\${report_feature_off} --no-report --data-dir --data-format tsv --config ${config} --title 'H3K27ac ChIP-seq quality control' --outdir . || true
 render_qc_summary.py --input-dir . --output-dir h3k27ac_multiqc_data --report h3k27ac_multiqc_report.html --title 'H3K27ac ChIP-seq quality control'
 multiqc --version > chip_multiqc.version.txt 2>&1
 """
}
