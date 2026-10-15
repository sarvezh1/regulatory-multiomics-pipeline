process BUILD_PROVENANCE {
    tag 'RNA provenance'
    cache false
    publishDir "${params.outdir}/rna/provenance", mode: 'copy', overwrite: true

    input:
    path samplesheet
    path reference_report
    path differential_status
    path differential_design
    path version_files, stageAs: 'versions??/*'
    path containers
    val nextflow_version
    val design_formula
    val contrast
    val run_de
    val execution_profile

    output:
    path 'rna.provenance.json', emit: manifest
    path 'rna.analysis_status.tsv', emit: statuses
    path 'provenance.version.txt', emit: version

    script:
    """
    build_provenance.py --samplesheet ${samplesheet} --reference ${reference_report} \
        --differential-status ${differential_status} --differential-design ${differential_design} \
        --versions ${version_files.join(' ')} --containers ${containers} \
        --nextflow-version '${nextflow_version}' --project-version '${params.project_version}' \
        --profile '${execution_profile}' --design '${design_formula}' \
        --contrast '${contrast}' --run-de '${run_de}' --reference-name '${params.reference_name}' \
        --reference-scope '${params.reference_scope}' --fasta-source-url '${params.fasta_source_url}' \
        --gtf-source-url '${params.gtf_source_url}' --manifest rna.provenance.json \
        --statuses rna.analysis_status.tsv
    python --version > provenance.version.txt 2>&1
    """
}
