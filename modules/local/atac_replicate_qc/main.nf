process ATAC_COUNT_REPLICATE_QC {
    tag 'ATAC count replicate QC'
    publishDir "${params.outdir}/atac/qc/replicates", mode: 'copy', overwrite: true
    input: path counts
    output: path 'atac.sample_correlation.tsv', emit: correlation; path 'atac.pca.tsv', emit: pca; path 'atac_replicate_qc.version.txt', emit: version
    script:
    """
    atac_replicate_qc.R --counts ${counts} --correlation atac.sample_correlation.tsv --pca atac.pca.tsv
    Rscript --version > atac_replicate_qc.version.txt 2>&1
    """
}

process ATAC_PEAK_JACCARD {
    tag 'ATAC peak Jaccard'
    publishDir "${params.outdir}/atac/qc/replicates", mode: 'copy', overwrite: true
    input: path peaks
    output: path 'atac.peak_jaccard.tsv', emit: jaccard; path 'atac_peak_jaccard.version.txt', emit: version
    script:
    """
    files=(${peaks.join(' ')})
    printf 'sample_a\tsample_b\tintersection_bases\tunion_bases\tjaccard\n' > atac.peak_jaccard.tsv
    for ((i=0; i<\${#files[@]}; i++)); do
      for ((j=i; j<\${#files[@]}; j++)); do
        a=\$(basename "\${files[\$i]}" | sed 's/_peaks.*//')
        b=\$(basename "\${files[\$j]}" | sed 's/_peaks.*//')
        read intersection union jaccard ignored <<< "\$(bedtools jaccard -a "\${files[\$i]}" -b "\${files[\$j]}" | tail -1)"
        if [ "\$union" = '0' ]; then jaccard='NA'; fi
        printf '%s\t%s\t%s\t%s\t%s\n' "\$a" "\$b" "\$intersection" "\$union" "\$jaccard" >> atac.peak_jaccard.tsv
      done
    done
    bedtools --version > atac_peak_jaccard.version.txt
    """
}
