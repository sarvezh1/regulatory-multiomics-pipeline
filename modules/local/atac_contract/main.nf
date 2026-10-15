process BUILD_ATAC_CONSENSUS {
    tag 'ATAC consensus regions'
    publishDir "${params.outdir}/atac/contract", mode: 'copy', overwrite: true
    input: path peaks
    output: path 'atac.regions.tsv', emit: regions; path 'atac.region_peak_provenance.tsv', emit: provenance; path 'atac_consensus.version.txt', emit: version
    script:
    """
    build_atac_consensus.py --peaks ${peaks.join(' ')} --regions atac.regions.tsv --provenance atac.region_peak_provenance.tsv
    python --version > atac_consensus.version.txt 2>&1
    """
}

process BUILD_ATAC_MATRIX {
    tag 'ATAC fragment count matrix'
    publishDir "${params.outdir}/atac/counts", mode: 'copy', overwrite: true
    input: path regions; path fragments
    output: path 'atac.raw_counts.tsv', emit: matrix; path 'atac_matrix.version.txt', emit: version
    script:
    """
    build_atac_matrix.py --regions ${regions} --fragments ${fragments.join(' ')} --output atac.raw_counts.tsv
    python --version > atac_matrix.version.txt 2>&1
    """
}

process ADAPT_PROCESSED_ATAC {
    tag 'processed ATAC region counts'
    publishDir "${params.outdir}/atac/counts", mode: 'copy', overwrite: true, pattern: 'atac.raw_counts.tsv'
    publishDir "${params.outdir}/atac/contract", mode: 'copy', overwrite: true, pattern: 'atac.regions.tsv'
    publishDir "${params.outdir}/atac/provenance", mode: 'copy', overwrite: true, pattern: 'atac.processed_input.provenance.json'
    input: path matrix; path metadata
    output: path 'atac.raw_counts.tsv', emit: counts; path 'atac.regions.tsv', emit: regions; path 'atac.processed_input.provenance.json', emit: provenance; path 'adapt_atac.version.txt', emit: version
    script:
    """
    adapt_atac_matrix.py --matrix ${matrix} --metadata ${metadata} --output atac.raw_counts.tsv --regions atac.regions.tsv --provenance atac.processed_input.provenance.json
    python --version > adapt_atac.version.txt 2>&1
    """
}
