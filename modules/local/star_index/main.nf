process STAR_INDEX {
    tag 'STAR index'
    publishDir "${params.outdir}/rna/provenance/reference_index", mode: 'copy', overwrite: true

    input:
    path fasta
    path gtf
    path validation_marker

    output:
    path 'star_index', emit: index
    path 'star.version.txt', emit: version

    script:
    """
    mkdir star_index
    STAR --runMode genomeGenerate \
        --runThreadN ${task.cpus} \
        --genomeDir star_index \
        --genomeFastaFiles ${fasta} \
        --sjdbGTFfile ${gtf} \
        --sjdbOverhang ${params.star_sjdb_overhang} \
        --genomeSAindexNbases ${params.star_sa_index_nbases} \
        --genomeChrBinNbits ${params.star_chr_bin_nbits}
    STAR --version > star.version.txt
    """
}
