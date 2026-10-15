process STAR_ALIGN {
    tag "${meta.library_id}"
    publishDir "${params.outdir}/rna/alignment", mode: 'copy', overwrite: true, pattern: '*.{bam,out,tab}'

    input:
    tuple val(meta), path(reads), path(fastp_json), path(fastp_html)
    path star_index

    output:
    tuple val(meta), path("${meta.library_id}.Aligned.sortedByCoord.out.bam"), path("${meta.library_id}.Log.final.out"), path("${meta.library_id}.Log.out"), path("${meta.library_id}.SJ.out.tab"), emit: alignment
    path 'star.version.txt', emit: version

    script:
    def readArgs = meta.library_layout == 'paired' ? "${reads[0]} ${reads[1]}" : "${reads[0]}"
    """
    STAR --runThreadN ${task.cpus} \
        --genomeDir ${star_index} \
        --readFilesIn ${readArgs} \
        --readFilesCommand zcat \
        --outTmpDir /tmp/star_work \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix ${meta.library_id}.
    STAR --version > star.version.txt
    """
}
