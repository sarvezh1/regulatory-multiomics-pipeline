process FEATURECOUNTS {
    tag "${meta.library_id}"
    publishDir "${params.outdir}/rna/counts/featurecounts", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(bam), path(log_final), path(log_out), path(sj)
    path gtf

    output:
    tuple val(meta), path("${meta.library_id}.featureCounts.tsv"), path("${meta.library_id}.featureCounts.tsv.summary"), emit: counts
    path 'featurecounts.version.txt', emit: version

    script:
    def strand = [unstranded: 0, forward: 1, reverse: 2][meta.strandedness]
    def paired = meta.library_layout == 'paired' ? '-p --countReadPairs' : ''
    """
    featureCounts -T ${task.cpus} -a ${gtf} -o ${meta.library_id}.featureCounts.tsv \
        -t exon -g gene_id -s ${strand} ${paired} ${bam}
    featureCounts -v > featurecounts.version.raw 2>&1
    awk 'NF { print; exit }' featurecounts.version.raw > featurecounts.version.txt
    """
}
