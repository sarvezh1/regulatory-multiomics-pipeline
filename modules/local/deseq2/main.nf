process RNA_DESEQ2 {
    tag 'RNA differential expression'
    publishDir "${params.outdir}/rna/differential", mode: 'copy', overwrite: true
    publishDir "${params.outdir}/rna/counts", mode: 'copy', overwrite: true, pattern: 'rna.normalized_counts.tsv'

    input:
    path matrix
    path metadata
    path annotations
    val design_formula
    val contrast
    val run_de

    output:
    path 'rna.normalized_counts.tsv', emit: normalized
    path 'rna.differential.complete.tsv', emit: differential
    path 'rna.differential.status.json', emit: status
    path 'rna.design.json', emit: design
    path 'deseq2.version.txt', emit: version

    script:
    """
    run_deseq2.R --counts ${matrix} --metadata ${metadata} --annotations ${annotations} \
        --design '${design_formula}' --contrast '${contrast}' --run-de '${run_de}' \
        --normalized rna.normalized_counts.tsv --results rna.differential.complete.tsv \
        --status rna.differential.status.json --design-output rna.design.json
    Rscript -e 'cat("R=", R.version.string, "\\nBioconductor=3.22\\nDESeq2=", as.character(packageVersion("DESeq2")), "\\n", sep="")' > deseq2.version.txt
    """
}
