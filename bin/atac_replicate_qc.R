#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
get_arg <- function(key) { i <- match(key,args); if(is.na(i) || i==length(args)) stop(paste("missing",key)); args[i+1] }
counts <- read.delim(get_arg("--counts"), check.names=FALSE, stringsAsFactors=FALSE)
sample_cols <- setdiff(names(counts), c("region_id","chromosome","start","end","coordinate_system"))
mat <- as.matrix(counts[,sample_cols,drop=FALSE]); storage.mode(mat) <- "numeric"
cor_mat <- if(nrow(mat) == 0) matrix(NA_real_,length(sample_cols),length(sample_cols),dimnames=list(sample_cols,sample_cols)) else if(length(sample_cols) > 1) cor(mat, method="pearson", use="pairwise.complete.obs") else matrix(1,1,1,dimnames=list(sample_cols,sample_cols))
cor_out <- data.frame(sample_id=rownames(cor_mat), cor_mat, check.names=FALSE)
write.table(cor_out, get_arg("--correlation"), sep="\t", quote=FALSE, row.names=FALSE, na="NA")
pca_out <- data.frame(sample_id=sample_cols, PC1=NA_real_, PC2=NA_real_, variance_PC1=NA_real_, variance_PC2=NA_real_)
if(length(sample_cols) >= 2) {
  keep <- apply(mat,1,var,na.rm=TRUE) > 0
  if(sum(keep) >= 2) {
    p <- prcomp(t(log2(mat[keep,,drop=FALSE]+1)), center=TRUE, scale.=FALSE)
    variance <- p$sdev^2/sum(p$sdev^2)
    pca_out$PC1 <- p$x[,1]
    if(ncol(p$x) >= 2) pca_out$PC2 <- p$x[,2]
    pca_out$variance_PC1 <- variance[1]
    if(length(variance) >= 2) pca_out$variance_PC2 <- variance[2]
  }
}
write.table(pca_out, get_arg("--pca"), sep="\t", quote=FALSE, row.names=FALSE, na="NA")
