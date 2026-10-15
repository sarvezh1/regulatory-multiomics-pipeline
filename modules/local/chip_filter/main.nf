process CHIP_FILTER {
 tag "${meta.library_id}"; publishDir "${params.outdir}/h3k27ac/qc/alignment", mode:'copy', overwrite:true, pattern:'*.{flagstat,stats,tsv}'
 input: tuple val(meta), path(sam), path(log)
 output: tuple val(meta), path("${meta.library_id}.preblacklist.bam"), path("${meta.library_id}.preblacklist.evidence.bed"), path("${meta.library_id}.filter_metrics.tsv"), emit:prepared; tuple val(meta), path("${meta.library_id}.prefilter.flagstat"), path("${meta.library_id}.prefilter.stats"), emit:qc; path 'chip_samtools.version.txt', emit:version
 script:
 def prep = meta.library_layout == 'paired' ? "samtools view -b ${sam} | samtools sort -n -o namesort.bam; samtools fixmate -m namesort.bam fixmate.bam; samtools sort -o positionsort.bam fixmate.bam" : "samtools view -b ${sam} | samtools sort -o positionsort.bam"
 def filt = meta.library_layout == 'paired' ? "samtools view -b -f 2 -F 3844 -q ${params.chip_mapq} marked.bam > ${meta.library_id}.preblacklist.bam" : "samtools view -b -F 3844 -q ${params.chip_mapq} marked.bam > ${meta.library_id}.preblacklist.bam"
 def evidence = meta.library_layout == 'paired' ? "samtools view ${meta.library_id}.preblacklist.bam | awk 'BEGIN{OFS=\"\\t\"} and(\$2,64) && \$9>0 {print \$3,\$4-1,\$4-1+\$9,\$1,60,\".\"}'" : "samtools view ${meta.library_id}.preblacklist.bam | awk 'BEGIN{OFS=\"\\t\"} {print \$3,\$4-1,\$4-1+length(\$10),\$1,60,and(\$2,16)?\"-\":\"+\"}'"
 def countflag = meta.library_layout == 'paired' ? '-f 64' : ''
 def proper = meta.library_layout == 'paired' ? "samtools view -c -f 66 -F 2308 marked.bam" : "samtools view -c -F 2308 marked.bam"
 def highq = meta.library_layout == 'paired' ? "samtools view -c -f 66 -F 2308 -q ${params.chip_mapq} marked.bam" : "samtools view -c -F 2308 -q ${params.chip_mapq} marked.bam"
 """
 ${prep}
 samtools markdup -s positionsort.bam marked.bam 2> ${meta.library_id}.markdup.stats
 samtools index marked.bam
 samtools flagstat marked.bam > ${meta.library_id}.prefilter.flagstat
 samtools stats marked.bam > ${meta.library_id}.prefilter.stats
 total=\$(samtools view -c ${countflag} marked.bam); mapped=\$(samtools view -c ${countflag} -F 4 marked.bam); proper=\$(${proper}); highq=\$(${highq}); dup=\$(samtools view -c ${countflag} -f 1024 marked.bam)
 if samtools view -H marked.bam | grep '^@SQ' | grep -Fq 'SN:${params.chip_mito_contig}'; then mt_available=1; mito=\$(samtools view -c ${countflag} -F 2308 -q ${params.chip_mapq} marked.bam '${params.chip_mito_contig}'); else mt_available=0; mito=0; fi
 samtools view marked.bam | awk 'BEGIN{OFS="\t"} ${meta.library_layout == 'paired' ? 'and($2,64)' : '1'} {k=\$3":"\$4":"${meta.library_layout == 'paired' ? '$9' : 'and($2,16)'}; n[k]++} END{for(k in n) print n[k],k}' > positions.txt
 unique=\$(wc -l < positions.txt); singleton=\$(awk '\$1==1{n++}END{print n+0}' positions.txt); doubleton=\$(awk '\$1==2{n++}END{print n+0}' positions.txt)
 nrf=\$(awk -v u=\$unique -v t=\$mapped 'BEGIN{if(t)printf "%.8f",u/t;else print "NA"}'); pbc1=\$(awk -v s=\$singleton -v u=\$unique 'BEGIN{if(u)printf "%.8f",s/u;else print "NA"}'); pbc2=\$(awk -v s=\$singleton -v d=\$doubleton 'BEGIN{if(d)printf "%.8f",s/d;else print "NA"}')
 ${filt}
 samtools view -b -e 'rname!="${params.chip_mito_contig}"' ${meta.library_id}.preblacklist.bam > nuclear.bam; mv nuclear.bam ${meta.library_id}.preblacklist.bam
 ${evidence} | sort -k1,1 -k2,2n > ${meta.library_id}.preblacklist.evidence.bed
 printf 'metric\tvalue\nall_reads_or_fragments\t%s\nmapped_reads_or_fragments\t%s\nproper_or_primary\t%s\nhigh_quality\t%s\nmitochondrial\t%s\nmitochondrial_reference_available\t%s\nduplicate\t%s\nunique_positions\t%s\nsingleton_positions\t%s\ndoubleton_positions\t%s\nNRF\t%s\nPBC1\t%s\nPBC2\t%s\n' "\$total" "\$mapped" "\$proper" "\$highq" "\$mito" "\$mt_available" "\$dup" "\$unique" "\$singleton" "\$doubleton" "\$nrf" "\$pbc1" "\$pbc2" > ${meta.library_id}.filter_metrics.tsv
 samtools --version | head -1 > chip_samtools.version.txt
 """
}

process CHIP_BLACKLIST_FILTER {
 tag "${meta.library_id}"
 input: tuple val(meta), path(bam), path(evidence), path(metrics, stageAs:'source.filter_metrics.tsv'); path blacklist
 output: tuple val(meta), path("${meta.library_id}.allowed_names.txt"), path("${meta.library_id}.evidence.bed"), path("${meta.library_id}.filter_metrics.tsv"), path(bam), emit:filtered; path 'chip_bedtools.version.txt', emit:version
 script:
 """
 bedtools intersect -v -a ${evidence} -b ${blacklist} > ${meta.library_id}.evidence.bed
 cut -f4 ${meta.library_id}.evidence.bed | sort -u > ${meta.library_id}.allowed_names.txt
 before=\$(wc -l < ${evidence}); after=\$(wc -l < ${meta.library_id}.evidence.bed)
 cp source.filter_metrics.tsv ${meta.library_id}.filter_metrics.tsv
 printf 'blacklist_input\t%s\nblacklist_retained\t%s\nblacklist_removed\t%s\n' "\$before" "\$after" "\$((before-after))" >> ${meta.library_id}.filter_metrics.tsv
 bedtools --version > chip_bedtools.version.txt
 """
}

process CHIP_FINAL_BAM {
 tag "${meta.library_id}"; publishDir "${params.outdir}/h3k27ac/alignment/libraries", mode:'copy', overwrite:true, pattern:'*.{bam,bai}'
 input: tuple val(meta), path(names), path(evidence), path(metrics), path(bam)
 output: tuple val(meta), path("${meta.library_id}.filtered.bam"), path("${meta.library_id}.filtered.bam.bai"), path(evidence), path(metrics), emit:filtered; path 'chip_final.version.txt', emit:version
 script:
 """
 if [ -s ${names} ]; then samtools view -N ${names} -b ${bam} | samtools sort -o ${meta.library_id}.filtered.bam; else samtools view -H ${bam} | samtools view -b -o ${meta.library_id}.filtered.bam; fi
 samtools index ${meta.library_id}.filtered.bam
 samtools --version | head -1 > chip_final.version.txt
 """
}
