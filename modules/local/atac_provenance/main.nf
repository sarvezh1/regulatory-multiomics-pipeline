process BUILD_ATAC_PROVENANCE {
    tag 'ATAC provenance'
    cache false
    publishDir "${params.outdir}/atac/provenance", mode: 'copy', overwrite: true
    input: path samplesheet; path validation_report; path reference_report; path blacklist; path tss; path processing_status; path da_status; path versions, stageAs: 'versions??/*'; path containers; val nextflow_version; val profile; val design; val contrast; val effective_genome_size
    output: path 'atac.provenance.json', emit: manifest; path 'atac.analysis_status.tsv', emit: status; path 'atac_provenance.version.txt', emit: version
    script:
    """
    cat ${versions.join(' ')} > atac.software_versions.txt
    build_atac_provenance.py --samplesheet ${samplesheet} --validation-report ${validation_report} --reference-report ${reference_report} --blacklist ${blacklist} --tss ${tss} --status ${processing_status} --da-status ${da_status} --versions atac.software_versions.txt --containers ${containers} --nextflow-version '${nextflow_version}' --profile '${profile}' --design '${design}' --contrast '${contrast}' --effective-genome-size '${effective_genome_size}' --output atac.provenance.json
    python -c 'import json; d=json.load(open("atac.provenance.json")); print("component\\tstatus"); [print(x["component"]+"\\t"+x["status"]) for x in d["statuses"]]' > atac.analysis_status.tsv
    python --version > atac_provenance.version.txt 2>&1
    """
}
