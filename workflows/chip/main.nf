include { VALIDATE_CHIP_SAMPLESHEET } from '../../modules/local/chip_validate'
include { VALIDATE_REFERENCE } from '../../modules/local/validate_reference'
include { CHIP_BOWTIE2_INDEX } from '../../modules/local/chip_reference'
include { CHIP_FASTQC_RAW; CHIP_FASTQC_TRIMMED } from '../../modules/local/chip_fastqc'
include { CHIP_FASTP } from '../../modules/local/chip_fastp'
include { CHIP_BOWTIE2_ALIGN } from '../../modules/local/chip_align'
include { CHIP_FILTER; CHIP_BLACKLIST_FILTER; CHIP_FINAL_BAM } from '../../modules/local/chip_filter'
include { CHIP_MERGE_SAMPLE } from '../../modules/local/chip_sample'
include { MACS3_H3K27AC } from '../../modules/local/chip_peaks'
include { CHIP_SIGNAL } from '../../modules/local/chip_signal'
include { CHIP_QC } from '../../modules/local/chip_qc'
include { BUILD_CHIP_CONSENSUS; BUILD_CHIP_MATRIX; ADAPT_PROCESSED_CHIP } from '../../modules/local/chip_contract'
include { CHIP_DIFFBIND } from '../../modules/local/chip_diff'
include { CHIP_PROCESSING_STATUS } from '../../modules/local/chip_status'
include { MULTIQC_CHIP } from '../../modules/local/chip_multiqc'
include { BUILD_CHIP_PROVENANCE } from '../../modules/local/chip_provenance'

workflow CHIP_WORKFLOW {
 take: samplesheet; fasta; gtf; blacklist; design; contrast; run_diff; route
 main:
 VALIDATE_CHIP_SAMPLESHEET(samplesheet)
 VALIDATE_REFERENCE(fasta,gtf,params.reference_scope)
 rows=VALIDATE_CHIP_SAMPLESHEET.out.validated.splitCsv(header:true,sep:'\t')
 if(route=='raw') {
   raw_reads=rows.filter{it.input_route=='raw'}.map{r->tuple(r as Map,r.read2?[file(r.read1),file(r.read2)]:[file(r.read1)])}
   CHIP_BOWTIE2_INDEX(fasta,VALIDATE_CHIP_SAMPLESHEET.out.raw_marker)
   CHIP_FASTQC_RAW(raw_reads)
   CHIP_FASTP(raw_reads)
   CHIP_FASTQC_TRIMMED(CHIP_FASTP.out.trimmed)
   CHIP_BOWTIE2_ALIGN(CHIP_FASTP.out.trimmed.map{m,r,j,h->tuple(m,r)},CHIP_BOWTIE2_INDEX.out.index)
   CHIP_FILTER(CHIP_BOWTIE2_ALIGN.out.alignment)
   CHIP_BLACKLIST_FILTER(CHIP_FILTER.out.prepared,blacklist)
   CHIP_FINAL_BAM(CHIP_BLACKLIST_FILTER.out.filtered)
   grouped=CHIP_FINAL_BAM.out.filtered.map{m,b,bi,e,met->tuple(m.sample_id,m,b,bi,e,met)}.groupTuple(by:0)
   CHIP_MERGE_SAMPLE(grouped)
   treat=CHIP_MERGE_SAMPLE.out.samples.filter{sid,m,b,bi,e,met->m[0].assay=='H3K27ac_ChIP'}.map{sid,m,b,bi,e,met->tuple(m[0].control_sample_id,sid,m,b,bi,e,met)}
   ctrl=CHIP_MERGE_SAMPLE.out.samples.filter{sid,m,b,bi,e,met->m[0].assay=='ChIP_input'}.map{sid,m,b,bi,e,met->tuple(sid,m,b,bi,e,met)}
   pairs=treat.join(ctrl,by:0).map{cid,sid,tm,tb,tbi,te,tmet,cm,cb,cbi,ce,cmet->tuple(sid,tm,tb,tbi,te,tmet,cid,cm,cb,cbi,ce,cmet)}
   MACS3_H3K27AC(pairs)
   CHIP_SIGNAL(MACS3_H3K27AC.out.peaks)
   CHIP_QC(MACS3_H3K27AC.out.peaks)
   peak_files=MACS3_H3K27AC.out.peaks.map{sid,tm,tb,tbi,te,tmet,cid,cm,cb,cbi,ce,cmet,p,s,x,pu,la,fe,l->p}.collect()
   evidence_files=MACS3_H3K27AC.out.peaks.map{sid,tm,tb,tbi,te,tmet,cid,cm,cb,cbi,ce,cmet,p,s,x,pu,la,fe,l->te}.collect()
   BUILD_CHIP_CONSENSUS(peak_files)
   BUILD_CHIP_MATRIX(BUILD_CHIP_CONSENSUS.out.regions,evidence_files)
   counts_ch=BUILD_CHIP_MATRIX.out.matrix
   CHIP_PROCESSING_STATUS(route,peak_files)
   qc_inputs=CHIP_FASTQC_RAW.out.reports.map{m,z,h->[z,h]}.mix(CHIP_FASTQC_TRIMMED.out.reports.map{m,z,h->[z,h]}).mix(CHIP_FASTP.out.trimmed.map{m,r,j,h->[j,h]}).mix(CHIP_BOWTIE2_ALIGN.out.alignment.map{m,s,l->[l]}).mix(CHIP_FILTER.out.qc.map{m,f,s->[f,s]}).flatten().collect()
   raw_versions=CHIP_BOWTIE2_INDEX.out.version.mix(CHIP_FASTQC_RAW.out.version).mix(CHIP_FASTP.out.version).mix(CHIP_FASTQC_TRIMMED.out.version).mix(CHIP_BOWTIE2_ALIGN.out.version).mix(CHIP_FILTER.out.version).mix(CHIP_BLACKLIST_FILTER.out.version).mix(CHIP_FINAL_BAM.out.version).mix(CHIP_MERGE_SAMPLE.out.version).mix(MACS3_H3K27AC.out.version).mix(CHIP_SIGNAL.out.version).mix(CHIP_QC.out.version).mix(BUILD_CHIP_CONSENSUS.out.version).mix(BUILD_CHIP_MATRIX.out.version)
 } else {
   matrix=rows.filter{it.input_route=='processed'}.map{r->file(r.processed_input)}.first()
   ADAPT_PROCESSED_CHIP(matrix,VALIDATE_CHIP_SAMPLESHEET.out.validated)
   counts_ch=ADAPT_PROCESSED_CHIP.out.counts
   CHIP_PROCESSING_STATUS(route,blacklist)
   qc_inputs=VALIDATE_CHIP_SAMPLESHEET.out.report.collect()
   raw_versions=ADAPT_PROCESSED_CHIP.out.version
 }
 CHIP_DIFFBIND(counts_ch,VALIDATE_CHIP_SAMPLESHEET.out.validated,design,contrast,run_diff)
 MULTIQC_CHIP(qc_inputs,file("${projectDir}/assets/multiqc_config.yml"))
 versions=VALIDATE_CHIP_SAMPLESHEET.out.version.mix(VALIDATE_REFERENCE.out.version).mix(raw_versions).mix(CHIP_DIFFBIND.out.version).mix(CHIP_PROCESSING_STATUS.out.version).mix(MULTIQC_CHIP.out.version).collect()
   BUILD_CHIP_PROVENANCE(VALIDATE_CHIP_SAMPLESHEET.out.validated,VALIDATE_CHIP_SAMPLESHEET.out.report,VALIDATE_REFERENCE.out.report,blacklist,CHIP_PROCESSING_STATUS.out.status,CHIP_DIFFBIND.out.status,versions,file("${projectDir}/assets/tool_containers.tsv"),workflow.nextflow.version,workflow.profile?:'standard',design,contrast,params.chip_effective_genome_size,params.chip_mapq,params.chip_se_validation_extsize)
}
