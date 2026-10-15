process ATAQV {
    tag "${sample_id}"
    publishDir "${params.outdir}/atac/qc/ataqv", mode: 'copy', overwrite: true
    input: tuple val(sample_id), val(metas), path(bam), path(bai), path(fragments), path(peaks), path(summits), path(xls), path(log), path(blacklist_metrics), path(filter_metrics); path tss
    output: tuple val(sample_id), path("${sample_id}.ataqv.json"), emit: metrics; path 'ataqv.version.txt', emit: version
    script:
    """
    ataqv --peak-file ${peaks} --tss-file ${tss} --metrics-file ${sample_id}.ataqv.json --name ${sample_id} human ${bam}
    (ataqv --version || true) > ataqv.version.txt 2>&1
    """
}

process ATAC_CUSTOM_QC {
    tag "${sample_id}"
    publishDir "${params.outdir}/atac/qc/metrics", mode: 'copy', overwrite: true
    input: tuple val(sample_id), val(metas), path(bam), path(bai), path(fragments), path(peaks), path(summits), path(xls), path(log), path(blacklist_metrics), path(filter_metrics); path tss
    output: tuple val(sample_id), path("${sample_id}.atac_qc.json"), path("${sample_id}.fragment_lengths.tsv"), emit: metrics; path 'atac_qc.version.txt', emit: version
    script:
    """
    atac_qc.py --fragments ${fragments} --pre-blacklist-metrics ${blacklist_metrics[0]} --filter-metrics ${filter_metrics[0]} --peaks ${peaks} --tss ${tss} --output ${sample_id}.atac_qc.json --lengths ${sample_id}.fragment_lengths.tsv
    python --version > atac_qc.version.txt 2>&1
    """
}
