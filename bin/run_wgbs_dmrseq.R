#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(dmrseq); library(bsseq); library(GenomicRanges)})
args <- commandArgs(trailingOnly=TRUE)
getarg <- function(name, default=NULL) { i<-match(name,args); if(is.na(i)) return(default); args[i+1] }
input <- getarg('--input'); metadata_file <- getarg('--metadata'); design <- getarg('--design'); contrast <- getarg('--contrast',''); run_dm <- tolower(getarg('--run-dm','false'))=='true'; min_cov <- as.integer(getarg('--min-coverage','5')); min_samples <- as.integer(getarg('--min-samples-per-group','2')); cutoff <- as.numeric(getarg('--cutoff','0.1')); max_perms <- as.integer(getarg('--max-perms','10')); output <- getarg('--output'); status_file <- getarg('--status'); design_file <- getarg('--design-output')
write_status <- function(status,reason,numerator='',denominator='',n_num=0,n_den=0,n_cpg=0,n_dmr=0){writeLines(jsonlite::toJSON(list(component='wgbs_differential',status=status,reason=reason,design_formula=design,contrast=contrast,numerator=numerator,denominator=denominator,independent_n_numerator=n_num,independent_n_denominator=n_den,minimum_coverage=min_cov,minimum_samples_per_group=min_samples,cpgs_analyzed=n_cpg,regions_returned=n_dmr,method='dmrseq',dmrseq_version=as.character(packageVersion('dmrseq')),bsseq_version=as.character(packageVersion('bsseq'))),auto_unbox=TRUE,pretty=TRUE),status_file)}
empty_output <- function(){writeLines('dmr_id\tchromosome\tstart\tend\tcoordinate_system\tn_cpg\tmean_methylation_numerator\tmean_methylation_denominator\tdelta_methylation\tdelta_methylation_pp\ttest_statistic\tpvalue\tqvalue\tcontrast\tmethod',output)}
empty_output()
writeLines('sample_id',design_file)
meta <- read.delim(metadata_file,check.names=FALSE,stringsAsFactors=FALSE)
long <- read.delim(input,check.names=FALSE,stringsAsFactors=FALSE,na.strings=c('NA',''))
if(!run_dm){write_status('NOT_RUN_DISABLED','Differential methylation disabled'); quit(save='no')}
if(contrast==''){write_status('NOT_RUN_NO_CONTRAST','No explicit contrast supplied'); quit(save='no')}
if(all(is.na(long$coverage))){write_status('NOT_AVAILABLE_INPUT','Count-preserving CpG coverage is unavailable'); quit(save='no')}
parts <- strsplit(contrast,',',fixed=TRUE)[[1]]; if(length(parts)!=3) stop('contrast must be variable,numerator,denominator')
variable<-parts[1]; numerator<-parts[2]; denominator<-parts[3]
terms <- trimws(strsplit(sub('^~','',design),'\\+')[[1]])
if(!all(terms %in% names(meta)) || !(variable %in% terms)) {write_status('NOT_RUN_INVALID_DESIGN','Design variables or contrast variable are invalid'); quit(save='no')}
meta <- meta[meta[[variable]] %in% c(numerator,denominator),,drop=FALSE]
if(!all(c(numerator,denominator) %in% meta[[variable]])){write_status('NOT_RUN_INVALID_DESIGN','Contrast levels absent',numerator,denominator); quit(save='no')}
n_num<-length(unique(meta$sample_id[meta[[variable]]==numerator])); n_den<-length(unique(meta$sample_id[meta[[variable]]==denominator]))
if(n_num<2 || n_den<2){write_status('NOT_RUN_INSUFFICIENT_REPLICATION','Each group requires at least two independent biological samples',numerator,denominator,n_num,n_den); quit(save='no')}
for(term in terms) meta[[term]]<-factor(meta[[term]])
meta[[variable]]<-relevel(meta[[variable]],ref=denominator)
mm<-model.matrix(as.formula(design),meta); if(qr(mm)$rank<ncol(mm)){write_status('NOT_RUN_INVALID_DESIGN','Design matrix is rank deficient',numerator,denominator,n_num,n_den); quit(save='no')}
write.table(cbind(sample_id=meta$sample_id,mm),design_file,sep='\t',quote=FALSE,row.names=FALSE)
long<-long[long$sample_id %in% meta$sample_id,,drop=FALSE]
if(!any(!is.na(long$coverage) & long$coverage>0)){write_status('COMPLETE_EMPTY_RESULT','No covered CpGs are available for the contrast',numerator,denominator,n_num,n_den,0,0); quit(save='no')}
key<-paste(long$chromosome,long$start,long$end,sep=':'); loci<-unique(long[,c('cpg_id','chromosome','start','end')]); loci$key<-paste(loci$chromosome,loci$start,loci$end,sep=':')
keep<-logical(nrow(loci)); names(keep)<-loci$key
for(k in loci$key){ z<-long[key==k & !is.na(long$coverage),]; keep[k]<-all(vapply(c(numerator,denominator),function(level) sum(z$sample_id %in% meta$sample_id[meta[[variable]]==level] & z$coverage>=min_cov)>=min_samples,logical(1))) }
loci<-loci[keep[loci$key],,drop=FALSE]
if(!nrow(loci)){write_status('COMPLETE_EMPTY_RESULT','No CpGs met the predeclared coverage rule',numerator,denominator,n_num,n_den,0,0); quit(save='no')}
sample_names<-meta$sample_id; M<-matrix(0,nrow=nrow(loci),ncol=length(sample_names),dimnames=list(loci$key,sample_names)); Cov<-M
for(i in seq_len(nrow(loci))) for(j in seq_along(sample_names)){ z<-long[key==loci$key[i] & long$sample_id==sample_names[j] & !is.na(long$coverage) & !is.na(long$methylated_count),]; if(nrow(z)){M[i,j]<-z$methylated_count[1]; Cov[i,j]<-z$coverage[1]} }
bs<-BSseq(chr=loci$chromosome,pos=loci$start+1,M=M,Cov=Cov,sampleNames=sample_names); pData(bs)<-S4Vectors::DataFrame(meta,row.names=meta$sample_id)
adjust<-setdiff(terms,variable); result<-dmrseq(bs=bs,testCovariate=variable,adjustCovariate=if(length(adjust)) adjust else NULL,cutoff=cutoff,minNumRegion=3,smooth=TRUE,minInSpan=3,maxPerms=max_perms,BPPARAM=BiocParallel::SerialParam(),verbose=FALSE)
if(length(result)==0){write_status('COMPLETE_EMPTY_RESULT','dmrseq completed without candidate regions',numerator,denominator,n_num,n_den,nrow(loci),0); quit(save='no')}
mc<-as.data.frame(mcols(result)); pick<-function(candidates,default=NA_real_){hit<-candidates[candidates %in% names(mc)]; if(length(hit)) mc[[hit[1]]] else rep(default,nrow(mc))}
rows<-vector('list',length(result))
for(i in seq_along(result)){inside<-loci$chromosome==as.character(seqnames(result)[i]) & loci$start>=start(result)[i]-1 & loci$end<=end(result)[i]; vals_num<-M[inside,meta[[variable]]==numerator,drop=FALSE]/Cov[inside,meta[[variable]]==numerator,drop=FALSE]; vals_den<-M[inside,meta[[variable]]==denominator,drop=FALSE]/Cov[inside,meta[[variable]]==denominator,drop=FALSE]; mn<-mean(vals_num,na.rm=TRUE); md<-mean(vals_den,na.rm=TRUE); delta<-mn-md; token<-paste(as.character(seqnames(result)[i]),start(result)[i]-1,end(result)[i],contrast,'dmrseq',sep=':'); id<-paste0('DMR_',substr(digest::digest(token,algo='sha256',serialize=FALSE),1,16)); rows[[i]]<-data.frame(dmr_id=id,chromosome=as.character(seqnames(result)[i]),start=start(result)[i]-1,end=end(result)[i],coordinate_system='0-based-half-open',n_cpg=sum(inside),mean_methylation_numerator=mn,mean_methylation_denominator=md,delta_methylation=delta,delta_methylation_pp=delta*100,test_statistic=pick(c('stat','testStat','area'))[i],pvalue=pick(c('pval','pvalue'))[i],qvalue=pick(c('qval','qvalue'))[i],contrast=contrast,method='dmrseq') }
write.table(do.call(rbind,rows),output,sep='\t',quote=FALSE,row.names=FALSE,na='NA'); write_status('COMPLETE','dmrseq completed',numerator,denominator,n_num,n_den,nrow(loci),length(result))
