process BISMARK_ALIGN {
    tag "${meta.library_id}"
    publishDir "${params.outdir}/wgbs/alignment", mode: 'copy', overwrite: true, pattern: '*report.txt'
    input: tuple val(meta), path(read1), path(read2), path(trim_report); path(index)
    output:
    tuple val(meta), path('*.bam'), path('*report.txt'), path(trim_report), emit: alignment
    path 'bismark_align.version.txt', emit: version
    script:
    def nonDirectional = meta.bisulfite_directionality == 'non_directional' ? '--non_directional' : ''
    """
    bismark align --bowtie2 --genome ${index} -1 ${read1} -2 ${read2} --output_dir . --basename ${meta.library_id} -p ${task.cpus} ${nonDirectional}
    { bismark --version; bowtie2 --version | head -1; } > bismark_align.version.txt 2>&1
    """
}
