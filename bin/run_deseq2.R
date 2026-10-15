#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
getarg <- function(name) { i <- match(paste0('--',name), args); if (is.na(i) || i==length(args)) stop(paste('missing',name)); args[i+1] }
counts_path<-getarg('counts'); metadata_path<-getarg('metadata'); annotations_path<-getarg('annotations')
design_text<-getarg('design'); contrast_text<-getarg('contrast'); run_de<-tolower(getarg('run-de'))=='true'
normalized_path<-getarg('normalized'); results_path<-getarg('results'); status_path<-getarg('status'); design_path<-getarg('design-output')
suppressPackageStartupMessages(library(DESeq2))
counts_in<-read.delim(counts_path,check.names=FALSE,stringsAsFactors=FALSE)
meta<-read.delim(metadata_path,check.names=FALSE,stringsAsFactors=FALSE)
ann<-read.delim(annotations_path,check.names=FALSE,stringsAsFactors=FALSE)
sample_names<-setdiff(colnames(counts_in),c('gene_id_original','gene_id_base','gene_name'))
if (!setequal(sample_names,meta$sample_id)) stop('count matrix and metadata samples differ')
meta<-meta[match(sample_names,meta$sample_id),,drop=FALSE]; rownames(meta)<-meta$sample_id
count_matrix<-as.matrix(counts_in[,sample_names,drop=FALSE]); storage.mode(count_matrix)<-'integer'; rownames(count_matrix)<-counts_in$gene_id_original
dds_norm<-DESeqDataSetFromMatrix(countData=count_matrix,colData=meta,design=~1)
if (nrow(dds_norm)>0 && any(rowSums(counts(dds_norm))>0)) { dds_norm<-estimateSizeFactors(dds_norm); norm<-counts(dds_norm,normalized=TRUE) } else norm<-count_matrix
norm_out<-data.frame(gene_id_original=counts_in$gene_id_original,gene_id_base=counts_in$gene_id_base,gene_name=counts_in$gene_name,norm,check.names=FALSE)
write.table(norm_out,normalized_path,sep='\t',quote=FALSE,row.names=FALSE,na='NA')
empty_results<-data.frame(gene_id_original=character(),gene_id_base=character(),gene_name=character(),baseMean=numeric(),log2FoldChange=numeric(),lfcSE=numeric(),stat=numeric(),pvalue=numeric(),padj=numeric())
status<-'COMPLETE'; reason<-'analysis completed'; compared_n<-list(); vars<-all.vars(as.formula(design_text)); contrast_parts<-if(nchar(contrast_text)) strsplit(contrast_text,',',fixed=TRUE)[[1]] else character()
if (!run_de) { status<-'NOT_RUN_DISABLED'; reason<-'differential analysis disabled explicitly'
} else if (length(contrast_parts)!=3) { status<-'NOT_RUN_NO_CONTRAST'; reason<-'no explicit contrast supplied'
} else if (!all(vars %in% colnames(meta))) { status<-'NOT_RUN_INVALID_DESIGN'; reason<-'design variable absent from sample metadata'
} else {
  variable<-contrast_parts[1]; test<-contrast_parts[2]; reference<-contrast_parts[3]
  if (!(variable %in% colnames(meta)) || !(variable %in% vars) || !all(c(test,reference) %in% unique(meta[[variable]]))) {
    status<-'NOT_RUN_INVALID_DESIGN'; reason<-'contrast variable is absent from the design or its requested levels are absent from sample metadata'
  } else {
    compared_n<-as.list(table(meta[[variable]][meta[[variable]] %in% c(test,reference)]))
    if (any(unlist(compared_n)<2)) { status<-'NOT_RUN_INSUFFICIENT_REPLICATION'; reason<-'each compared group requires at least two independent biological samples'
    } else {
      for (v in vars) meta[[v]]<-factor(meta[[v]])
      mm<-try(model.matrix(as.formula(design_text),data=meta),silent=TRUE)
      if (inherits(mm,'try-error') || qr(mm)$rank<ncol(mm)) { status<-'NOT_RUN_INVALID_DESIGN'; reason<-'design matrix is not full rank'
      } else {
        dds<-DESeqDataSetFromMatrix(countData=count_matrix,colData=meta,design=as.formula(design_text))
        dds<-DESeq(dds,quiet=TRUE)
        res<-results(dds,contrast=c(variable,test,reference),alpha=0.05)
        z<-as.data.frame(res); z$gene_id_original<-rownames(z); z<-merge(ann,z,by='gene_id_original',all.y=TRUE,sort=FALSE)
        z<-z[,c('gene_id_original','gene_id_base','gene_name','baseMean','log2FoldChange','lfcSE','stat','pvalue','padj')]
        write.table(z,results_path,sep='\t',quote=FALSE,row.names=FALSE,na='NA')
        if (sum(!is.na(z$padj) & z$padj<0.05)==0) { status<-'COMPLETE_EMPTY_RESULT'; reason<-'valid analysis completed with zero genes below the declared FDR threshold' }
      }
    }
  }
}
if (!file.exists(results_path)) write.table(empty_results,results_path,sep='\t',quote=FALSE,row.names=FALSE,na='NA')
esc<-function(x) gsub('"','\\\\"',gsub('\\\\','\\\\\\\\',x))
status_json<-sprintf('{\n  "status": "%s",\n  "reason": "%s",\n  "alpha": 0.05,\n  "method": "DESeq2",\n  "independent_biological_samples": %d\n}\n',status,esc(reason),nrow(meta))
writeLines(status_json,status_path)
levels_text<-paste(sprintf('"%s": [%s]',vars,sapply(vars,function(v) paste(sprintf('"%s"',sort(unique(meta[[v]]))),collapse=', '))),collapse=', ')
design_json<-sprintf('{\n  "design_formula": "%s",\n  "contrast": "%s",\n  "samples_included": [%s],\n  "samples_excluded": [],\n  "factor_levels": {%s},\n  "filtering_rule": "DESeq2 independent filtering",\n  "alpha": 0.05\n}\n',esc(design_text),esc(contrast_text),paste(sprintf('"%s"',meta$sample_id),collapse=', '),levels_text)
writeLines(design_json,design_path)
