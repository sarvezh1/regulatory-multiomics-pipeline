#!/usr/bin/env bash
set -euo pipefail
out_dir=${1:-tests/integration/reference}
mkdir -p "$out_dir"
tmp_dir=$(mktemp -d)
trap 'rm -rf "$tmp_dir"' EXIT
fasta_url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000022.11&rettype=fasta&retmode=text'
gtf_url='https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.annotation.gtf.gz'
blacklist_url='https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz'
curl -L --fail --retry 3 "$fasta_url" -o "$tmp_dir/NC_000022.11.fa"
curl -L --fail --retry 3 "$gtf_url" -o "$tmp_dir/gencode.v50.annotation.gtf.gz"
curl -L --fail --retry 3 "$blacklist_url" -o "$out_dir/ENCFF356LFX.GRCh38.blacklist.bed.gz"
awk 'NR==1{sub(/^>.*/,">chr22")} {print}' "$tmp_dir/NC_000022.11.fa" > "$out_dir/GRCh38.p14.primary_assembly.chr22.fa"
gzip -cd "$tmp_dir/gencode.v50.annotation.gtf.gz" | awk 'BEGIN{FS=OFS="\t"} /^#/ || $1=="chr22"' > "$out_dir/gencode.v50.annotation.chr22.gtf"
gzip -cd "$out_dir/ENCFF356LFX.GRCh38.blacklist.bed.gz" | awk 'BEGIN{FS=OFS="\t"} $1=="chr22"' > "$out_dir/ENCFF356LFX.GRCh38.blacklist.chr22.bed"
{
  printf 'artifact\tsha256\tsize_bytes\n'
  printf 'NC_000022.11.source.fa\t%s\t%s\n' "$(sha256sum "$tmp_dir/NC_000022.11.fa" | cut -d' ' -f1)" "$(stat -c%s "$tmp_dir/NC_000022.11.fa")"
  printf 'gencode.v50.annotation.gtf.gz\t%s\t%s\n' "$(sha256sum "$tmp_dir/gencode.v50.annotation.gtf.gz" | cut -d' ' -f1)" "$(stat -c%s "$tmp_dir/gencode.v50.annotation.gtf.gz")"
  for f in "$out_dir/GRCh38.p14.primary_assembly.chr22.fa" "$out_dir/gencode.v50.annotation.chr22.gtf"; do
    printf '%s\t%s\t%s\n' "$(basename "$f")" "$(sha256sum "$f" | cut -d' ' -f1)" "$(stat -c%s "$f")"
  done
  for f in "$out_dir/ENCFF356LFX.GRCh38.blacklist.bed.gz" "$out_dir/ENCFF356LFX.GRCh38.blacklist.chr22.bed"; do
    printf '%s\t%s\t%s\n' "$(basename "$f")" "$(sha256sum "$f" | cut -d' ' -f1)" "$(stat -c%s "$f")"
  done
} > "$out_dir/checksums.tsv"
