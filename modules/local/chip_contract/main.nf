process BUILD_CHIP_CONSENSUS {
 tag 'H3K27ac consensus'; publishDir "${params.outdir}/h3k27ac/contract", mode:'copy', overwrite:true
 input: path peaks
 output: path 'h3k27ac.regions.tsv', emit:regions; path 'h3k27ac.region_peak_provenance.tsv', emit:provenance; path 'chip_consensus.version.txt', emit:version
 script:
 """
 build_chip_consensus.py --peaks ${peaks.join(' ')} --regions h3k27ac.regions.tsv --provenance h3k27ac.region_peak_provenance.tsv
 python --version > chip_consensus.version.txt 2>&1
 """
}
process BUILD_CHIP_MATRIX {
 tag 'H3K27ac count matrix'; publishDir "${params.outdir}/h3k27ac/counts", mode:'copy', overwrite:true
 input: path regions; path evidence
 output: path 'h3k27ac.raw_counts.tsv', emit:matrix; path 'chip_matrix.version.txt', emit:version
 script:
 """
 build_chip_matrix.py --regions ${regions} --evidence ${evidence.join(' ')} --output h3k27ac.raw_counts.tsv
 python --version > chip_matrix.version.txt 2>&1
 """
}
process ADAPT_PROCESSED_CHIP {
 tag 'processed H3K27ac counts'; publishDir "${params.outdir}/h3k27ac/counts", mode:'copy', overwrite:true, pattern:'h3k27ac.raw_counts.tsv'; publishDir "${params.outdir}/h3k27ac/contract", mode:'copy', overwrite:true, pattern:'h3k27ac.regions.tsv'; publishDir "${params.outdir}/h3k27ac/provenance", mode:'copy', overwrite:true, pattern:'h3k27ac.processed_input.provenance.json'
 input: path matrix; path metadata
 output: path 'h3k27ac.raw_counts.tsv', emit:counts; path 'h3k27ac.regions.tsv', emit:regions; path 'h3k27ac.processed_input.provenance.json', emit:provenance; path 'adapt_chip.version.txt', emit:version
 script:
 """
 adapt_processed_chip.py --matrix ${matrix} --metadata ${metadata} --counts h3k27ac.raw_counts.tsv --regions h3k27ac.regions.tsv --provenance h3k27ac.processed_input.provenance.json
 python --version > adapt_chip.version.txt 2>&1
 """
}
