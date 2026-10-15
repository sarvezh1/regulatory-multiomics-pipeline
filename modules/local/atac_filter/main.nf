process ATAC_SAMTOOLS_PREP {
    tag "${meta.library_id}"
    publishDir "${params.outdir}/atac/qc/alignment", mode: 'copy', overwrite: true, pattern: '*.{flagstat,stats,tsv}'
    input: tuple val(meta), path(sam), path(bowtie_log)
    output: tuple val(meta), path("${meta.library_id}.preblacklist.bam"), path("${meta.library_id}.preblacklist.fragments.bed"), path("${meta.library_id}.filter_metrics.tsv"), emit: prepared; tuple val(meta), path("${meta.library_id}.prefilter.flagstat"), path("${meta.library_id}.prefilter.stats"), emit: qc; path 'samtools_atac.version.txt', emit: version
    script:
    """
    samtools view -@ ${task.cpus} -b ${sam} | samtools sort -@ ${task.cpus} -n -o namesort.bam
    samtools fixmate -@ ${task.cpus} -m namesort.bam fixmate.bam
    samtools sort -@ ${task.cpus} -o positionsort.bam fixmate.bam
    samtools markdup -@ ${task.cpus} -s positionsort.bam marked.bam 2> ${meta.library_id}.markdup.stats
    samtools index marked.bam
    samtools flagstat marked.bam > ${meta.library_id}.prefilter.flagstat
    samtools stats marked.bam > ${meta.library_id}.prefilter.stats
    total=\$(samtools view -c -f 64 marked.bam)
    mapped=\$(samtools view -c -f 64 -F 4 marked.bam)
    proper=\$(samtools view -c -f 66 -F 2308 marked.bam)
    highq=\$(samtools view -c -f 66 -F 2308 -q ${params.atac_mapq} marked.bam)
    if samtools view -H marked.bam | grep '^@SQ' | grep -Fq 'SN:${params.atac_mito_contig}'; then mt_available=1; mito=\$(samtools view -c -f 66 -F 2308 -q ${params.atac_mapq} marked.bam '${params.atac_mito_contig}'); else mt_available=0; mito=0; fi
    dup=\$(samtools view -c -f 1088 marked.bam)
    samtools view marked.bam | awk 'BEGIN{OFS="\t"} and(\$2,64) && !and(\$2,4) && \$9>0 {print \$3,\$4-1,\$4-1+\$9}' | sort | uniq -c > fragment_positions.txt
    unique=\$(wc -l < fragment_positions.txt)
    singleton=\$(awk '\$1==1{n++} END{print n+0}' fragment_positions.txt)
    doubleton=\$(awk '\$1==2{n++} END{print n+0}' fragment_positions.txt)
    nrf=\$(awk -v u="\$unique" -v t="\$mapped" 'BEGIN{if(t) printf "%.8f",u/t; else print "NA"}')
    pbc1=\$(awk -v s="\$singleton" -v u="\$unique" 'BEGIN{if(u) printf "%.8f",s/u; else print "NA"}')
    pbc2=\$(awk -v s="\$singleton" -v d="\$doubleton" 'BEGIN{if(d) printf "%.8f",s/d; else print "NA"}')
    samtools view -@ ${task.cpus} -b -f 2 -F 3844 -q ${params.atac_mapq} marked.bam | samtools view -@ ${task.cpus} -b -e 'rname!="${params.atac_mito_contig}"' -o ${meta.library_id}.preblacklist.bam
    samtools index ${meta.library_id}.preblacklist.bam
    samtools sort -n -o final.namesort.bam ${meta.library_id}.preblacklist.bam
    samtools view final.namesort.bam | awk 'BEGIN{OFS="\t"} \$9>0 {print \$3,\$4-1,\$4-1+\$9,\$1,60,"."}' > ${meta.library_id}.preblacklist.fragments.bed
    printf 'metric\tvalue\nall_first_mates\t%s\nmapped_first_mates\t%s\nproper_primary\t%s\nhigh_quality\t%s\nmitochondrial\t%s\nmitochondrial_reference_available\t%s\nduplicate\t%s\nunique_fragment_positions\t%s\nsingleton_positions\t%s\ndoubleton_positions\t%s\nNRF\t%s\nPBC1\t%s\nPBC2\t%s\n' "\$total" "\$mapped" "\$proper" "\$highq" "\$mito" "\$mt_available" "\$dup" "\$unique" "\$singleton" "\$doubleton" "\$nrf" "\$pbc1" "\$pbc2" > ${meta.library_id}.filter_metrics.tsv
    samtools --version | head -1 > samtools_atac.version.txt
    """
}

process ATAC_BLACKLIST_FILTER {
    tag "${meta.library_id}"
    input: tuple val(meta), path(bam), path(fragments), path(filter_metrics); path blacklist
    output: tuple val(meta), path("${meta.library_id}.allowed_names.txt"), path("${meta.library_id}.fragments.bed"), path("${meta.library_id}.blacklist_metrics.json"), path(bam), path(filter_metrics), emit: filtered_names; path 'blacklist_filter.version.txt', emit: version
    script:
    """
    filter_atac_fragments.py --fragments ${fragments} --blacklist ${blacklist} --output ${meta.library_id}.fragments.bed --names ${meta.library_id}.allowed_names.txt --metrics ${meta.library_id}.blacklist_metrics.json
    python --version > blacklist_filter.version.txt 2>&1
    """
}

process ATAC_FINAL_BAM {
    tag "${meta.library_id}"
    publishDir "${params.outdir}/atac/alignment/libraries", mode: 'copy', overwrite: true, pattern: '*.{bam,bai}'
    input: tuple val(meta), path(names), path(fragments), path(blacklist_metrics), path(bam), path(filter_metrics)
    output: tuple val(meta), path("${meta.library_id}.filtered.bam"), path("${meta.library_id}.filtered.bam.bai"), path(fragments), path(blacklist_metrics), path(filter_metrics), emit: filtered; path 'samtools_final.version.txt', emit: version
    script:
    """
    if [ -s ${names} ]; then samtools view -@ ${task.cpus} -N ${names} -b ${bam} | samtools sort -@ ${task.cpus} -o ${meta.library_id}.filtered.bam; else samtools view -H ${bam} | samtools view -b -o ${meta.library_id}.filtered.bam; fi
    samtools index ${meta.library_id}.filtered.bam
    samtools --version | head -1 > samtools_final.version.txt
    """
}
