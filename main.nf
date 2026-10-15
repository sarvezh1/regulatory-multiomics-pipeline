nextflow.enable.dsl = 2

include { RNA_WORKFLOW } from './workflows/rna'
include { ATAC_WORKFLOW } from './workflows/atac'
include { WGBS_WORKFLOW } from './workflows/wgbs'
include { CHIP_WORKFLOW } from './workflows/chip'
include { INTEGRATION_WORKFLOW } from './workflows/integration'

workflow {
    if (params.help.toString().toBoolean() || params.pipeline_help.toString().toBoolean()) {
        log.info file("${projectDir}/assets/help.txt").text
    } else {
      def designPattern = ~/^~\s*[A-Za-z][A-Za-z0-9_]*(\s*\+\s*[A-Za-z][A-Za-z0-9_]*)*$/
      def contrastPattern = ~/^[A-Za-z][A-Za-z0-9_]*,[A-Za-z0-9_.-]+,[A-Za-z0-9_.-]+$/
      if (params.assay != 'INTEGRATION' && !params.input) error "--input is required"
      if (!params.fasta) error "--fasta is required"
      if (!params.gtf) error "--gtf is required"
      if (!(params.design_formula ==~ designPattern)) error "--design_formula must contain only '~ variable' or '~ covariate + variable' terms"
      if (params.contrast && !(params.contrast ==~ contrastPattern)) error "--contrast must be variable,test_level,reference_level"
      if (params.assay == 'RNA') {
        RNA_WORKFLOW(file(params.input, checkIfExists: true), file(params.fasta, checkIfExists: true), file(params.gtf, checkIfExists: true), params.design_formula, params.contrast ?: '', params.run_de.toString().toBoolean())
    } else if (params.assay == 'ATAC') {
        if (!params.blacklist) error "--blacklist is required for ATAC"
        if (!(params.atac_input_route in ['raw','processed'])) error "--atac_input_route must be raw or processed"
        ATAC_WORKFLOW(file(params.input, checkIfExists: true), file(params.fasta, checkIfExists: true), file(params.gtf, checkIfExists: true), file(params.blacklist, checkIfExists: true), params.design_formula, params.contrast ?: '', params.run_da.toString().toBoolean(), params.atac_input_route)
    } else if (params.assay == 'WGBS') {
        if (!(params.wgbs_input_route in ['raw','processed'])) error "--wgbs_input_route must be raw or processed"
        WGBS_WORKFLOW(file(params.input, checkIfExists: true), file(params.fasta, checkIfExists: true), file(params.gtf, checkIfExists: true), params.design_formula, params.contrast ?: '', params.run_dm.toString().toBoolean(), params.wgbs_input_route)
    } else if (params.assay == 'H3K27ac') {
        if (!params.blacklist) error "--blacklist is required for H3K27ac"
        if (!(params.chip_input_route in ['raw','processed'])) error "--chip_input_route must be raw or processed"
        CHIP_WORKFLOW(file(params.input, checkIfExists: true), file(params.fasta, checkIfExists: true), file(params.gtf, checkIfExists: true), file(params.blacklist, checkIfExists: true), params.design_formula, params.contrast ?: '', params.run_chip_diff.toString().toBoolean(), params.chip_input_route)
    } else if (params.assay == 'INTEGRATION') {
        if (!params.integration_manifest) error "--integration_manifest is required for INTEGRATION"
        if (!(params.integration_mode in ['baseline','differential','both'])) error "--integration_mode must be baseline, differential, or both"
        if (!(params.integration_match_type in ['EXACT_SAMPLE_MATCHED','CONDITION_MATCHED','REFERENCE_COMPARABLE','UNMATCHED'])) error "invalid --integration_match_type"
        INTEGRATION_WORKFLOW(file(params.integration_manifest, checkIfExists: true), file(params.fasta, checkIfExists: true), file(params.gtf, checkIfExists: true), params.integration_mode, params.integration_match_type, params.integration_fdr, params.distal_max_distance)
      } else {
        error "--assay must be exactly RNA, ATAC, WGBS, H3K27ac, or INTEGRATION; unsupported assays are rejected"
      }
    }
}
