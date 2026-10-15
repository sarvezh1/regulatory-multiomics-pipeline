process WGBS_DMRSEQ {
    tag 'WGBS differential methylation'
    publishDir "${params.outdir}/wgbs/differential", mode: 'copy', overwrite: true
    input: path long_table; path metadata; val design; val contrast; val run_dm
    output:
    path 'wgbs.dmrseq.complete.tsv', emit: results
    path 'wgbs.differential_status.json', emit: status
    path 'wgbs.design_matrix.tsv', emit: design
    path 'wgbs_dmrseq.version.txt', emit: version
    script:
    """
    run_wgbs_dmrseq.R --input ${long_table} --metadata ${metadata} --design '${design}' --contrast '${contrast}' --run-dm ${run_dm} --min-coverage ${params.wgbs_min_coverage} --min-samples-per-group ${params.wgbs_min_samples_per_group} --cutoff ${params.dmrseq_cutoff} --max-perms ${params.dmrseq_max_perms} --output wgbs.dmrseq.complete.tsv --status wgbs.differential_status.json --design-output wgbs.design_matrix.tsv
    Rscript -e "cat('R ',R.version[['major']],'.',R.version[['minor']],'; Bioconductor ',as.character(BiocManager::version()),'; dmrseq ',as.character(packageVersion('dmrseq')),'; bsseq ',as.character(packageVersion('bsseq')),'\\n',sep='')" > wgbs_dmrseq.version.txt
    """
}
