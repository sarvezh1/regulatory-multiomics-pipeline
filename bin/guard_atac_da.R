#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE); opt<-list(); i<-1
while(i<=length(args)){opt[[sub('^--','',args[i])]]<-args[i+1]; i<-i+2}
counts<-read.delim(opt$counts,check.names=FALSE,stringsAsFactors=FALSE); meta<-read.delim(opt$metadata,check.names=FALSE,stringsAsFactors=FALSE)
run_da<-tolower(opt$`run-da`)=='true'; design<-opt$design; contrast<-opt$contrast
status<-'READY'; reason<-'design guard passed'; vars<-trimws(strsplit(sub('^~','',design),'\\+')[[1]])
parts<-if(nzchar(contrast)) strsplit(contrast,',',fixed=TRUE)[[1]] else character()
samples<-names(counts)[-(1:4)]; meta<-meta[match(samples,meta$sample_id),,drop=FALSE]
if(!run_da){status<-'NOT_RUN_DISABLED';reason<-'differential accessibility disabled explicitly'} else if(length(parts)!=3){status<-'NOT_RUN_NO_CONTRAST';reason<-'run_da requested without an explicit contrast'} else if(!parts[1]%in%vars || !all(vars%in%names(meta)) || any(is.na(meta$sample_id)) || anyDuplicated(meta$sample_id)){status<-'NOT_RUN_INVALID_DESIGN';reason<-'invalid design fields or biological-sample metadata'} else {
  groups<-table(meta[[parts[1]]]); if(!all(parts[2:3]%in%names(groups)) || any(groups[parts[2:3]]<2)){status<-'NOT_RUN_INSUFFICIENT_REPLICATION';reason<-'each contrast group requires at least two independent biological samples'} else { for(v in vars) meta[[v]]<-factor(meta[[v]]); mm<-model.matrix(as.formula(design),meta); if(qr(mm)$rank<ncol(mm)){status<-'NOT_RUN_INVALID_DESIGN';reason<-'design matrix is rank deficient'} }
}
if(status=='READY' && nrow(counts)==0){status<-'COMPLETE_EMPTY_RESULT';reason<-'valid design has no consensus regions to test'}
if(status=='READY'){writeLines('ready',opt$ready);quit(status=0)}
missing_stat<-rep(NA_real_,nrow(counts))
result<-data.frame(region_id=counts$region_id,chromosome=counts$chromosome,start=counts$start,end=counts$end,baseMean=missing_stat,log2FoldChange=missing_stat,lfcSE=missing_stat,stat=missing_stat,pvalue=missing_stat,padj=missing_stat)
write.table(result,opt$results,sep='\t',quote=FALSE,row.names=FALSE,na='NA'); write.table(counts,opt$normalized,sep='\t',quote=FALSE,row.names=FALSE,na='NA')
esc<-function(x) gsub('"','\\\\"',x,fixed=TRUE)
writeLines(c('{',sprintf('  "status": "%s",',status),sprintf('  "reason": "%s",',esc(reason)),sprintf('  "design_formula": "%s",',esc(design)),sprintf('  "contrast": "%s",',esc(contrast)),sprintf('  "independent_biological_samples": %d',nrow(meta)),'}'),opt$status)
