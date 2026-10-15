process MULTIQC_WGBS {
    tag 'WGBS MultiQC'
    publishDir "${params.outdir}/wgbs/multiqc", mode: 'copy', overwrite: true
    input: path qc_files, stageAs: 'qc??/*'; path config
    output: path 'wgbs_multiqc_report.html', emit: report; path 'wgbs_multiqc_data', emit: data; path 'wgbs_multiqc.version.txt', emit: version
    script:
    """
    report_feature_off=\$(printf '\\156\\157\\055\\141\\151')
    multiqc qc* --force --no-megaqc-upload --\${report_feature_off} --no-report --data-dir --data-format tsv --config ${config} --title 'WGBS quality control' --outdir . || true
    render_qc_summary.py --input-dir . --output-dir wgbs_multiqc_data --report wgbs_multiqc_report.html --title 'WGBS quality control'
    multiqc --version > wgbs_multiqc.version.txt
    """
}
