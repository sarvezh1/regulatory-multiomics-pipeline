process ATAC_DIFFBIND {
    tag 'ATAC differential accessibility'
    publishDir "${params.outdir}/atac/differential", mode: 'copy', overwrite: true
    publishDir "${params.outdir}/atac/counts", mode: 'copy', overwrite: true, pattern: 'atac.normalized_counts.tsv'
    input: path counts; path metadata; val design; val contrast; val run_da
    output: path 'atac.differential.complete.tsv', emit: differential; path 'atac.normalized_counts.tsv', emit: normalized; path 'atac.differential.status.json', emit: status; path 'atac_da.version.txt', emit: version
    script:
    """
    guard_atac_da.R --counts ${counts} --metadata ${metadata} --design '${design}' --contrast '${contrast}' --run-da '${run_da}' --results atac.differential.complete.tsv --normalized atac.normalized_counts.tsv --status atac.differential.status.json --ready atac_da.ready
    if [ -s atac_da.ready ]; then
      run_atac_da.R --counts ${counts} --metadata ${metadata} --design '${design}' --contrast '${contrast}' --run-da true --results atac.differential.complete.tsv --normalized atac.normalized_counts.tsv --status atac.differential.status.json
    else
      Rscript --version > atac_da.version.txt 2>&1
      printf 'Bioconductor=3.22\nDiffBind=not_loaded_guarded\nDESeq2=not_loaded_guarded\n' >> atac_da.version.txt
    fi
    """
}
