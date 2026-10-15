process CHIP_DIFFBIND {
 tag 'H3K27ac differential enrichment'; publishDir "${params.outdir}/h3k27ac/differential", mode:'copy', overwrite:true
 input: path counts; path metadata; val design; val contrast; val run_diff
 output: path 'h3k27ac.differential.complete.tsv', emit:results; path 'h3k27ac.differential.status.json', emit:status; path 'chip_diff.version.txt', emit:version
 script:
 """
 run_chip_diff.R --counts ${counts} --metadata ${metadata} --design '${design}' --contrast '${contrast}' --run '${run_diff}' --output h3k27ac.differential.complete.tsv --status h3k27ac.differential.status.json
 """
}
