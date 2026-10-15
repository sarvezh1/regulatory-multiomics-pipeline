process WGBS_LIBRARY_CONTRACT {
    tag "${meta.library_id}"
    publishDir "${params.outdir}/wgbs/counts/libraries", mode: 'copy', overwrite: true, pattern: '*.tsv'
    input: tuple val(meta), path(coverage), path(mbias), path(split_report), path(chg), path(chh), path(dedup_report), path(alignment_report), path(trim_report); path(fasta)
    output:
    tuple val(meta), path("${meta.library_id}.cpg.tsv"), path("${meta.library_id}.nonreference_cpg.tsv"), path(mbias), path(split_report), path(chg), path(chh), path(dedup_report), path(alignment_report), path(trim_report), emit: contract
    path 'wgbs_library_contract.version.txt', emit: version
    script:
    """
    wgbs_cpg_contract.py --coverage ${coverage} --fasta ${fasta} --sample-id ${meta.sample_id} --library-id ${meta.library_id} --output ${meta.library_id}.cpg.tsv --rejected-output ${meta.library_id}.nonreference_cpg.tsv
    python --version > wgbs_library_contract.version.txt 2>&1
    """
}

process ADAPT_PROCESSED_WGBS {
    tag "${meta.library_id}"
    publishDir "${params.outdir}/wgbs/counts/libraries", mode: 'copy', overwrite: true, pattern: '*.cpg.tsv'
    input: tuple val(meta), path(processed); path(fasta)
    output: tuple val(meta), path("${meta.library_id}.cpg.tsv"), emit: contract
    path 'adapt_processed_wgbs.version.txt', emit: version
    script:
    """
    adapt_processed_wgbs.py --input ${processed} --fasta ${fasta} --sample-id ${meta.sample_id} --library-id ${meta.library_id} --type ${meta.processed_input_type} --output ${meta.library_id}.cpg.tsv
    python --version > adapt_processed_wgbs.version.txt 2>&1
    """
}

process AGGREGATE_WGBS_SAMPLE {
    tag "${sample_id}"
    publishDir "${params.outdir}/wgbs/counts/samples", mode: 'copy', overwrite: true, pattern: '*.cpg.tsv'
    input: tuple val(sample_id), path(contracts)
    output: tuple val(sample_id), path("${sample_id}.cpg.tsv"), emit: samples
    path 'aggregate_wgbs_sample.version.txt', emit: version
    script:
    """
    aggregate_wgbs_contracts.py --sample-id ${sample_id} --inputs ${contracts.join(' ')} --output ${sample_id}.cpg.tsv
    python --version > aggregate_wgbs_sample.version.txt 2>&1
    """
}

process BUILD_WGBS_UNIVERSE {
    tag 'WGBS CpG universe'
    publishDir "${params.outdir}/wgbs/counts", mode: 'copy', overwrite: true
    input: path contracts; path metadata
    output:
    path 'wgbs.cpg.long.tsv', emit: long
    path 'wgbs.coverage_summary.tsv', emit: summary
    path 'build_wgbs_universe.version.txt', emit: version
    script:
    """
    build_wgbs_universe.py --inputs ${contracts.join(' ')} --metadata ${metadata} --min-coverage ${params.wgbs_min_coverage} --output wgbs.cpg.long.tsv --summary wgbs.coverage_summary.tsv
    python --version > build_wgbs_universe.version.txt 2>&1
    """
}
