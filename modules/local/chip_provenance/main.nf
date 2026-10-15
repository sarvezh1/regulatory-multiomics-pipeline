process BUILD_CHIP_PROVENANCE {
 tag 'H3K27ac provenance'; cache false; publishDir "${params.outdir}/h3k27ac/provenance", mode:'copy', overwrite:true
 input: path samplesheet; path validation; path reference; path blacklist; path processing_status; path diff_status; path versions, stageAs:'versions??/*'; path containers; val nextflow_version; val profile; val design; val contrast; val effective_genome_size; val mapq; val se_validation_extsize
 output: path 'h3k27ac.provenance.json', emit:manifest; path 'h3k27ac.analysis_status.tsv', emit:status; path 'chip_provenance.version.txt', emit:version
 script:
 """
 cat ${versions.join(' ')} > chip.software_versions.txt
 build_chip_provenance.py --samplesheet ${samplesheet} --validation ${validation} --reference ${reference} --blacklist ${blacklist} --processing-status ${processing_status} --diff-status ${diff_status} --versions chip.software_versions.txt --containers ${containers} --nextflow-version '${nextflow_version}' --profile '${profile}' --design '${design}' --contrast '${contrast}' --effective-genome-size '${effective_genome_size}' --mapq '${mapq}' --se-validation-extsize '${se_validation_extsize}' --output h3k27ac.provenance.json
 python -c 'import json;d=json.load(open("h3k27ac.provenance.json"));print("component\\tstatus");[print(x["component"]+"\\t"+x["status"]) for x in d["statuses"]]' > h3k27ac.analysis_status.tsv
 python --version > chip_provenance.version.txt 2>&1
 """
}
