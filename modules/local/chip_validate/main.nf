process VALIDATE_CHIP_SAMPLESHEET {
 tag 'H3K27ac samplesheet'; publishDir "${params.outdir}/h3k27ac/provenance", mode:'copy', overwrite:true
 input: path samplesheet
 output: path 'chip.validated_samplesheet.tsv', emit:validated; path 'chip.samplesheet_validation.json', emit:report; path 'raw_inputs.present', optional:true, emit:raw_marker; path 'processed_inputs.present', optional:true, emit:processed_marker; path 'chip_validate.version.txt', emit:version
 script:
 """
 validate_chip_samplesheet.py --input ${samplesheet} --output chip.validated_samplesheet.tsv --report chip.samplesheet_validation.json
 python --version > chip_validate.version.txt 2>&1
 """
}
