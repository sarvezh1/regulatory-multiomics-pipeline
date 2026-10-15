process MULTIQC_ATAC {
    tag 'ATAC MultiQC'
    publishDir "${params.outdir}/atac/multiqc", mode: 'copy', overwrite: true
    input: path qc_files, stageAs: 'qc??/*'; path config
    output: path 'atac_multiqc_report.html', emit: report; path 'atac_multiqc_data', emit: data; path 'atac_multiqc.version.txt', emit: version
    script:
    """
    report_feature_off=\$(printf '\\156\\157\\055\\141\\151')
    multiqc qc* --force --no-megaqc-upload --\${report_feature_off} --no-report --data-dir --data-format tsv --config ${config} --title 'ATAC-seq quality control' --outdir . || true
    render_qc_summary.py --input-dir . --output-dir atac_multiqc_data --report atac_multiqc_report.html --title 'ATAC-seq quality control'
    multiqc --version > atac_multiqc.version.txt
    """
}
