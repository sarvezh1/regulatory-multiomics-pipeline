process RUN_REGULATORY_INTEGRATION {
    tag "${integration_mode} regulatory integration"
    publishDir params.outdir, mode: 'copy', overwrite: true

    input:
    path manifest
    path contract_inputs, stageAs: 'inputs/*'
    path reference_dir
    val integration_mode
    val integration_match_type
    val integration_fdr
    val distal_max_distance

    output:
    path 'integration', emit: results

    script:
    """
    stage_integration_manifest.py --manifest ${manifest} --input-dir inputs --output staged.integration_manifest.tsv
    integrate_regulatory.py \
      --manifest staged.integration_manifest.tsv \
      --promoters ${reference_dir}/integration.promoters.m2000_p500.tsv \
      --tss ${reference_dir}/integration.tss.tsv \
      --sensitivity-promoters ${reference_dir}/integration.promoters.m1000_p500.tsv ${reference_dir}/integration.promoters.m2000_p1000.tsv \
      --outdir integration \
      --mode '${integration_mode}' \
      --match-type '${integration_match_type}' \
      --fdr ${integration_fdr} \
      --distal-max ${distal_max_distance}
    cp -r ${reference_dir} integration/reference
    python --version > integration/provenance/integration_runtime.version.txt 2>&1
    """
}
