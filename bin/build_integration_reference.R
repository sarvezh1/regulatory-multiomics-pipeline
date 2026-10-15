#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
  library(rtracklayer)
  library(Biostrings)
})

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 4) stop("usage: build_integration_reference.R GTF FASTA OUTDIR ANNOTATION")
gtf_path <- args[[1]]; fasta_path <- args[[2]]; outdir <- args[[3]]; annotation <- args[[4]]
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

allowed <- c(paste0("chr", 1:22), "chrX", "chrY")
gtf <- import(gtf_path)
tx <- gtf[gtf$type == "transcript" & as.character(seqnames(gtf)) %in% allowed]
if (!length(tx)) stop("GTF contains no transcript features on primary nuclear chromosomes")
if (!all(c("gene_id", "transcript_id") %in% colnames(mcols(tx)))) stop("GTF transcript features require gene_id and transcript_id")

seqs <- readDNAStringSet(fasta_path)
names(seqs) <- sub(" .*", "", names(seqs))
seq_lengths <- width(seqs); names(seq_lengths) <- names(seqs)
missing_seq <- setdiff(unique(as.character(seqnames(tx))), names(seq_lengths))
if (length(missing_seq)) stop(paste("GTF chromosomes absent from FASTA:", paste(missing_seq, collapse=",")))

base_id <- function(x) sub("\\.[0-9]+$", "", x)
gene_original <- as.character(tx$gene_id)
gene_base <- base_id(gene_original)
gene_name <- if ("gene_name" %in% colnames(mcols(tx))) as.character(tx$gene_name) else rep("NA", length(tx))
strand_v <- as.character(strand(tx))
tss_one <- ifelse(strand_v == "+", start(tx), end(tx))
key <- paste(as.character(seqnames(tx)), tss_one, strand_v, gene_base, sep="\t")
groups <- split(seq_along(tx), key)

rows <- lapply(groups, function(ix) {
  chr <- as.character(seqnames(tx))[ix[1]]; pos1 <- tss_one[ix[1]]; st <- strand_v[ix[1]]
  gb <- gene_base[ix[1]]
  tids <- paste(sort(unique(as.character(tx$transcript_id[ix]))), collapse=";")
  originals <- paste(sort(unique(gene_original[ix])), collapse=";")
  gnames <- sort(unique(gene_name[ix][!is.na(gene_name[ix]) & gene_name[ix] != ""]))
  gname <- if (length(gnames)) paste(gnames, collapse=";") else "NA"
  tss_id <- paste("TSS", gb, chr, pos1, ifelse(st == "+", "PLUS", "MINUS"), sep="_")
  data.frame(tss_id=tss_id, gene_id_original=originals, gene_id_base=gb, gene_name=gname,
             transcript_ids=tids, chromosome=chr, tss_start=pos1-1L, tss_end=pos1,
             strand=st, coordinate_system="0-based-half-open", source_annotation=annotation,
             stringsAsFactors=FALSE)
})
tss <- do.call(rbind, rows)
ord <- order(factor(tss$chromosome, levels=allowed), tss$tss_start, tss$strand, tss$gene_id_base)
tss <- tss[ord, , drop=FALSE]

write_tsv <- function(x, path) write.table(x, path, sep="\t", quote=FALSE, row.names=FALSE, na="NA")
write_tsv(tss, file.path(outdir, "integration.tss.tsv"))

make_promoters <- function(upstream, downstream, label) {
  gr <- GRanges(seqnames=tss$chromosome,
                ranges=IRanges(start=tss$tss_start+1L, end=tss$tss_end),
                strand=tss$strand)
  seqinfo(gr) <- Seqinfo(seqnames=names(seq_lengths), seqlengths=as.integer(seq_lengths))
  pr <- trim(promoters(gr, upstream=upstream, downstream=downstream))
  out <- tss[, c("tss_id","gene_id_original","gene_id_base","gene_name","chromosome","strand","source_annotation")]
  out$promoter_start <- start(pr)-1L
  out$promoter_end <- end(pr)
  out$coordinate_system <- "0-based-half-open"
  out$upstream_bp <- upstream
  out$downstream_bp <- downstream
  out$promoter_id <- paste("PROM", out$tss_id, paste0("U", upstream), paste0("D", downstream), sep="_")
  out <- out[, c("promoter_id","tss_id","gene_id_original","gene_id_base","gene_name","chromosome",
                 "promoter_start","promoter_end","strand","coordinate_system","upstream_bp","downstream_bp","source_annotation")]
  write_tsv(out, file.path(outdir, paste0("integration.promoters.", label, ".tsv")))
}
make_promoters(2000L, 500L, "m2000_p500")
make_promoters(1000L, 500L, "m1000_p500")
make_promoters(2000L, 1000L, "m2000_p1000")

versions <- c(
  paste("R", as.character(getRversion())),
  paste("Bioconductor", "3.22"),
  paste("GenomicRanges", as.character(packageVersion("GenomicRanges"))),
  paste("IRanges", as.character(packageVersion("IRanges"))),
  paste("S4Vectors", as.character(packageVersion("S4Vectors"))),
  paste("rtracklayer", as.character(packageVersion("rtracklayer"))),
  paste("Biostrings", as.character(packageVersion("Biostrings")))
)
writeLines(versions, file.path(outdir, "integration_reference.versions.txt"))
