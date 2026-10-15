#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(DiffBind);library(DESeq2)})
args<-commandArgs(trailingOnly=TRUE); opt<-list(); i<-1
while(i<=length(args)){opt[[sub('^--','',args[i])]]<-args[i+1];i<-i+2}
counts<-read.delim(opt$counts,check.names=FALSE,stringsAsFactors=FALSE); meta<-read.delim(opt$metadata,check.names=FALSE,stringsAsFactors=FALSE)
meta<-meta[meta$assay=='H3K27ac_ChIP',,drop=FALSE]; run<-tolower(opt$run)=='true'; design<-opt$design; contrast<-opt$contrast
status<-'NOT_RUN_DISABLED';reason<-'differential H3K27ac enrichment disabled explicitly';vars<-trimws(strsplit(sub('^~','',design),'\\+')[[1]]); result<-data.frame()
if(run&&!nzchar(contrast)){status<-'NOT_RUN_NO_CONTRAST';reason<-'run requested without explicit contrast'} else if(run){
 p<-strsplit(contrast,',',fixed=TRUE)[[1]]; samples<-intersect(names(counts)[-(1:4)],meta$sample_id); meta<-meta[match(samples,meta$sample_id),,drop=FALSE]
 if(length(p)!=3||!p[1]%in%vars||!all(vars%in%names(meta))){status<-'NOT_RUN_INVALID_DESIGN';reason<-'invalid design or contrast'} else if(!all(c(p[2],p[3])%in%meta[[p[1]]])||any(table(meta[[p[1]]])[c(p[2],p[3])]<2)){status<-'NOT_RUN_INSUFFICIENT_REPLICATION';reason<-'each group requires at least two independent biological samples'} else {
  for(v in vars)meta[[v]]<-factor(meta[[v]]); mm<-model.matrix(as.formula(design),meta)
  if(qr(mm)$rank<ncol(mm)){status<-'NOT_RUN_INVALID_DESIGN';reason<-'design matrix is rank deficient'} else {
   mat<-as.matrix(counts[,samples,drop=FALSE]);storage.mode(mat)<-'integer';rownames(mat)<-counts$region_id; dba_obj<-NULL
   for(j in seq_along(samples)){ps<-data.frame(Chr=counts$chromosome,Start=counts$start,End=counts$end,Score=mat[,j]);dba_obj<-dba.peakset(dba_obj,ps,sampID=samples[j],condition=as.character(meta[[p[1]]][j]),replicate=as.integer(meta$biological_replicate[j]),peak.caller='counts',scoreCol=4)}
   if(nrow(dba.show(dba_obj))!=length(samples))stop('DiffBind sample framework lost a biological sample')
   dds<-DESeqDataSetFromMatrix(mat,meta,as.formula(design)); dds<-DESeq(dds,fitType='mean',quiet=TRUE); rr<-results(dds,contrast=c(p[1],p[2],p[3]),alpha=.05)
   result<-cbind(counts[,1:4],as.data.frame(rr)[counts$region_id,]);status<-'COMPLETE';reason<-paste0('DiffBind ',packageVersion('DiffBind'),' sample/region framework; DESeq2 ',packageVersion('DESeq2'),' count-model backend; positive log2FoldChange means greater H3K27ac enrichment in ',p[2],' than ',p[3])
  }
 }
}
if(status!='COMPLETE'){z<-rep(NA_real_,nrow(counts));result<-data.frame(counts[,1:4],baseMean=z,log2FoldChange=z,lfcSE=z,stat=z,pvalue=z,padj=z)}
write.table(result,opt$output,sep='\t',quote=FALSE,row.names=FALSE,na='NA');jsonlite::write_json(list(status=status,reason=reason,design_formula=design,contrast=contrast,independent_biological_samples=nrow(meta),DiffBind=as.character(packageVersion('DiffBind')),DESeq2=as.character(packageVersion('DESeq2')),control_handling='matched controls required for raw peak definition; no manual subtraction; processed route retains control metadata/counts'),opt$status,auto_unbox=TRUE,pretty=TRUE,na='string')
writeLines(c(paste0('R=',R.version.string),'Bioconductor=3.22',paste0('DiffBind=',packageVersion('DiffBind')),paste0('DESeq2=',packageVersion('DESeq2'))),'chip_diff.version.txt')
