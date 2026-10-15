process CHIP_BOWTIE2_INDEX {
 tag 'H3K27ac Bowtie2 index'; publishDir "${params.outdir}/h3k27ac/provenance/reference", mode:'copy', overwrite:true
 input: path fasta; path marker
 output: path 'bowtie2_index', emit:index; path 'chip_bowtie2_build.version.txt', emit:version
 script:
 """
 mkdir bowtie2_index
 bowtie2-build --threads ${task.cpus} ${fasta} bowtie2_index/reference
 bowtie2-build --version | head -1 > chip_bowtie2_build.version.txt
 """
}
