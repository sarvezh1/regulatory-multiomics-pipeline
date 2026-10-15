process FEATURECOUNTS_TO_CONTRACT {
    tag "${meta.library_id}"
    publishDir "${params.outdir}/rna/counts/per_library", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(counts), path(summary)
    path gtf

    output:
    tuple val(meta), path("${meta.library_id}.rna_counts.tsv"), emit: contract
    path 'featurecounts_contract.version.txt', emit: version

    script:
    """
    featurecounts_to_contract.py --counts ${counts} --gtf ${gtf} \
        --sample-id '${meta.sample_id}' --library-id '${meta.library_id}' \
        --source-genome-build '${params.reference_name}' --source-annotation '${params.annotation_name}' \
        --source-processing-pipeline 'workflow_raw_rna' \
        --source-file-accession 'NOT_APPLICABLE' --source-file-checksum 'NOT_APPLICABLE' \
        --coordinate-normalization-status 'not_required' --gene-id-mapping-status 'not_required' \
        --output ${meta.library_id}.rna_counts.tsv
    python --version > featurecounts_contract.version.txt 2>&1
    """
}

process ADAPT_PROCESSED_COUNTS {
    tag "${meta.library_id}"
    publishDir "${params.outdir}/rna/counts/per_library", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(counts)

    output:
    tuple val(meta), path("${meta.library_id}.rna_counts.tsv"), emit: contract
    path 'processed_adapter.version.txt', emit: version

    script:
    """
    adapt_processed_counts.py --input ${counts} --input-type '${meta.processed_input_type}' \
        --sample-id '${meta.sample_id}' --library-id '${meta.library_id}' \
        --source-genome-build '${meta.source_genome_build}' --source-annotation '${meta.source_annotation}' \
        --source-processing-pipeline '${meta.source_processing_pipeline}' \
        --source-file-accession '${meta.source_file_accession}' --source-file-checksum '${meta.source_file_checksum}' \
        --coordinate-normalization-status '${meta.coordinate_normalization_status}' \
        --gene-id-mapping-status '${meta.gene_id_mapping_status}' \
        --output ${meta.library_id}.rna_counts.tsv
    python --version > processed_adapter.version.txt 2>&1
    """
}

process MERGE_RNA_COUNTS {
    tag 'RNA count matrix'
    publishDir "${params.outdir}/rna/counts", mode: 'copy', overwrite: true

    input:
    path contracts
    path samplesheet

    output:
    path 'rna.raw_counts.tsv', emit: matrix
    path 'rna.sample_metadata.tsv', emit: metadata
    path 'rna.gene_annotations.tsv', emit: annotations
    path 'rna.technical_replicate_manifest.tsv', emit: technical_manifest
    path 'merge_contract.version.txt', emit: version

    script:
    """
    merge_rna_counts.py --samplesheet ${samplesheet} --contracts ${contracts.join(' ')} \
        --matrix rna.raw_counts.tsv --metadata rna.sample_metadata.tsv \
        --annotations rna.gene_annotations.tsv --technical-manifest rna.technical_replicate_manifest.tsv
    python --version > merge_contract.version.txt 2>&1
    """
}
