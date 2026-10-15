#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(DiffBind); library(DESeq2)})
args <- commandArgs(trailingOnly=TRUE); opt <- list(); i<-1
while(i<=length(args)){opt[[sub('^--','',args[i])]]<-args[i+1]; i<-i+2}
counts <- read.delim(opt$counts,check.names=FALSE,stringsAsFactors=FALSE)
meta <- read.delim(opt$metadata,check.names=FALSE,stringsAsFactors=FALSE)
run_da <- tolower(opt$`run-da`)=='true'; design_text<-opt$design; contrast_text<-opt$contrast
status<-'NOT_RUN_DISABLED'; reason<-'differential accessibility disabled explicitly'; result<-data.frame()
vars<-trimws(strsplit(sub('^~','',design_text),'\\+')[[1]])
valid_name<-function(x) grepl('^[A-Za-z][A-Za-z0-9_]*$',x)
if(run_da && !nzchar(contrast_text)){status<-'NOT_RUN_NO_CONTRAST'; reason<-'run_da requested without an explicit contrast'} else if(run_da){
  parts<-strsplit(contrast_text,',',fixed=TRUE)[[1]]
  if(length(parts)!=3 || !all(vapply(c(vars,parts[1]),valid_name,logical(1))) || !parts[1] %in% vars || !all(vars %in% names(meta))){status<-'NOT_RUN_INVALID_DESIGN'; reason<-'invalid design or contrast fields'} else {
    samples<-names(counts)[-(1:4)]; meta<-meta[match(samples,meta$sample_id),,drop=FALSE]
    if(any(is.na(meta$sample_id)) || anyDuplicated(meta$sample_id)){status<-'NOT_RUN_INVALID_DESIGN'; reason<-'metadata does not provide one row per count column'} else {
      test<-parts[2]; ref<-parts[3]; groups<-table(meta[[parts[1]]])
      if(!all(c(test,ref)%in%names(groups)) || any(groups[c(test,ref)]<2)){status<-'NOT_RUN_INSUFFICIENT_REPLICATION'; reason<-'each contrast group requires at least two independent biological samples'} else {
        for(v in vars) meta[[v]]<-factor(meta[[v]])
        mm<-model.matrix(as.formula(design_text),meta)
        if(qr(mm)$rank<ncol(mm)){status<-'NOT_RUN_INVALID_DESIGN'; reason<-'design matrix is rank deficient'} else {
          mat<-as.matrix(counts[,samples,drop=FALSE]); storage.mode(mat)<-'integer'; rownames(mat)<-counts$region_id
          dba_object<-NULL
          for(j in seq_along(samples)) {
            peakset<-data.frame(Chr=counts$chromosome,Start=counts$start,End=counts$end,Score=mat[,j])
            replicate_value<-if('biological_replicate' %in% names(meta)) as.integer(meta$biological_replicate[j]) else j
            dba_object<-dba.peakset(dba_object,peakset,sampID=samples[j],condition=as.character(meta[[parts[1]]][j]),replicate=replicate_value,peak.caller='counts',scoreCol=4,bLowerScoreBetter=FALSE)
          }
          if(nrow(dba.show(dba_object)) != length(samples)) stop('DiffBind sample model did not preserve all biological samples')
          dds<-DESeqDataSetFromMatrix(mat,meta,as.formula(design_text))
          dispersion_method<-'mean'
          dds<-tryCatch(DESeq(dds,fitType='mean',quiet=TRUE),error=function(e){
            if(!grepl('gene-wise dispersion estimates',conditionMessage(e),fixed=TRUE)) stop(e)
            dispersion_method<<-'gene_wise_fallback'
            x<-estimateSizeFactors(dds); x<-estimateDispersionsGeneEst(x,quiet=TRUE); dispersions(x)<-mcols(x)$dispGeneEst; nbinomWaldTest(x,quiet=TRUE)
          })
          rr<-results(dds,contrast=c(parts[1],test,ref),alpha=0.05)
          result<-cbind(counts[,1:4],as.data.frame(rr)[counts$region_id,]); names(result)[names(result)=='lfcSE']<-'lfcSE'
          status<-'COMPLETE'; reason<-paste0('DiffBind ',packageVersion('DiffBind'),' sample/peak framework with DESeq2 ',packageVersion('DESeq2'),' common-region count-model backend; dispersion=',dispersion_method,'; positive log2FoldChange means greater accessibility in ',test,' than ',ref)
          norm<-counts[,1:4]; norm<-cbind(norm,as.data.frame(counts(dds,normalized=TRUE))[counts$region_id,,drop=FALSE]); write.table(norm,opt$normalized,sep='\t',quote=FALSE,row.names=FALSE,na='NA')
        }
      }
    }
  }
}
if(status!='COMPLETE'){
  missing_stat<-rep(NA_real_,nrow(counts))
  result<-data.frame(region_id=counts$region_id,chromosome=counts$chromosome,start=counts$start,end=counts$end,baseMean=missing_stat,log2FoldChange=missing_stat,lfcSE=missing_stat,stat=missing_stat,pvalue=missing_stat,padj=missing_stat)
  write.table(counts,opt$normalized,sep='\t',quote=FALSE,row.names=FALSE,na='NA')
}
write.table(result,opt$results,sep='\t',quote=FALSE,row.names=FALSE,na='NA')
jsonlite::write_json(list(status=status,reason=reason,design_formula=design_text,contrast=contrast_text,independent_biological_samples=nrow(meta),DiffBind=as.character(packageVersion('DiffBind')),DESeq2=as.character(packageVersion('DESeq2'))),opt$status,auto_unbox=TRUE,pretty=TRUE,na='string')
writeLines(c(paste0('R=',R.version.string),'Bioconductor=3.22',paste0('DiffBind=',packageVersion('DiffBind')),paste0('DESeq2=',packageVersion('DESeq2'))),'atac_da.version.txt')
