process BUILD_WGBS_PROVENANCE {
    tag 'WGBS provenance'
    cache false
    publishDir "${params.outdir}/wgbs/provenance", mode: 'copy', overwrite: true
    input: path samplesheet; path validation_report; path reference_report; path reference_checksums; path processing_status; path dm_status; path versions, stageAs: 'versions??/*'; path containers; val nextflow_version; val profile; val design; val contrast
    output: path 'wgbs.provenance.json', emit: manifest; path 'wgbs.analysis_status.tsv', emit: status; path 'wgbs_provenance.version.txt', emit: version
    script:
    """
    cat ${versions.join(' ')} > wgbs.software_versions.txt
    build_wgbs_provenance.py --samplesheet ${samplesheet} --validation-report ${validation_report} --reference-report ${reference_report} --reference-checksums ${reference_checksums} --processing-status ${processing_status} --dm-status ${dm_status} --versions wgbs.software_versions.txt --containers ${containers} --nextflow-version '${nextflow_version}' --profile '${profile}' --design '${design}' --contrast '${contrast}' --min-coverage ${params.wgbs_min_coverage} --ignore ${params.wgbs_ignore} --ignore-r2 ${params.wgbs_ignore_r2} --ignore-3prime ${params.wgbs_ignore_3prime} --ignore-3prime-r2 ${params.wgbs_ignore_3prime_r2} --output wgbs.provenance.json
    python -c 'import json; d=json.load(open("wgbs.provenance.json")); print("component\\tstatus"); [print(x["component"]+"\\t"+x["status"]) for x in d["statuses"]]' > wgbs.analysis_status.tsv
    python --version > wgbs_provenance.version.txt 2>&1
    """
}
