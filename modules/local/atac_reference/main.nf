process BOWTIE2_INDEX {
    tag 'Bowtie2 index'
    publishDir "${params.outdir}/atac/provenance/reference", mode: 'copy', overwrite: true
    input: path fasta; path raw_marker
    output: path 'bowtie2_index', emit: index; path 'bowtie2_build.version.txt', emit: version
    script:
    """
    mkdir bowtie2_index
    bowtie2-build --threads ${task.cpus} ${fasta} bowtie2_index/reference
    bowtie2-build --version | head -1 > bowtie2_build.version.txt
    """
}

process DERIVE_ATAC_TSS {
    tag 'deduplicated TSS resource'
    publishDir "${params.outdir}/atac/provenance/reference", mode: 'copy', overwrite: true
    input: path gtf; path raw_marker
    output: path 'atac.tss.bed', emit: bed; path 'atac.tss.provenance.tsv', emit: provenance; path 'derive_tss.version.txt', emit: version
    script:
    """
    derive_atac_tss.py --gtf ${gtf} --bed atac.tss.bed --provenance atac.tss.provenance.tsv
    python --version > derive_tss.version.txt 2>&1
    """
}
