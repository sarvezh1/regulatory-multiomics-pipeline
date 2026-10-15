process BISMARK_REFERENCE {
    tag 'Bismark reference'
    publishDir "${params.outdir}/wgbs/provenance/reference_index", mode: 'copy', overwrite: true
    input: path fasta; path raw_marker
    output:
    path 'bismark_reference', emit: index
    path 'bismark_reference.checksums.tsv', emit: checksums
    path 'bismark_reference.version.txt', emit: version
    script:
    """
    mkdir bismark_reference
    cp ${fasta} bismark_reference/reference.fa
    bismark prepare --bowtie2 --parallel ${task.cpus} bismark_reference
    printf 'file\tsha256\n' > bismark_reference.checksums.tsv
    find bismark_reference -type f -print0 | sort -z | while IFS= read -r -d '' f; do printf '%s\t%s\n' "\${f#bismark_reference/}" "\$(sha256sum "\$f" | cut -d' ' -f1)"; done >> bismark_reference.checksums.tsv
    { bismark --version; bowtie2 --version | head -1; } > bismark_reference.version.txt 2>&1
    """
}

process WGBS_REFERENCE_CHECKSUM {
    tag 'WGBS reference checksum'
    input: path fasta
    output: path 'wgbs_reference.checksums.tsv', emit: checksums; path 'wgbs_reference_checksum.version.txt', emit: version
    script:
    """
    printf 'file\tsha256\n%s\t%s\n' '${fasta.name}' "\$(sha256sum ${fasta} | cut -d' ' -f1)" > wgbs_reference.checksums.tsv
    sha256sum --version | head -1 > wgbs_reference_checksum.version.txt
    """
}
