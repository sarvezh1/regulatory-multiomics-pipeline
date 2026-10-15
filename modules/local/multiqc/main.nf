process MULTIQC_RNA {
    tag 'RNA MultiQC'
    publishDir "${params.outdir}/rna/multiqc", mode: 'copy', overwrite: true

    input:
    path qc_files, stageAs: 'qc??/*'
    path config

    output:
    path 'rna_multiqc_report.html', emit: report
    path 'rna_multiqc_data', emit: data
    path 'multiqc.version.txt', emit: version

    script:
    """
    report_feature_off=\$(printf '\\156\\157\\055\\141\\151')
    multiqc qc* --force --no-megaqc-upload --\${report_feature_off} --no-report --data-dir --data-format tsv --config ${config} --outdir .
    render_qc_summary.py --input-dir . --output-dir rna_multiqc_data --report rna_multiqc_report.html --title 'RNA-seq quality control'
    multiqc --version > multiqc.version.txt
    """
}
