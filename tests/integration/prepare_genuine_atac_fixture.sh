#!/usr/bin/env bash
set -euo pipefail
out_dir=${1:-tests/integration/genuine_atac}
spots=${2:-100000}
mkdir -p "$out_dir"
for run in SRR13036031 SRR13036032; do
  source_url="https://sra-pub-run-odp.s3.amazonaws.com/sra/${run}/${run}"
  rm -f "$out_dir/${run}_1.fastq" "$out_dir/${run}_2.fastq" "$out_dir/${run}_1.fastq.gz" "$out_dir/${run}_2.fastq.gz"
  docker run --rm -v "$(pwd)/$out_dir:/data" -w /data \
    ncbi/sra-tools:3.4.1@sha256:4f048ba6754246e1ed7799a2c8ed469a38427f392e403ca08b93c4255104c9ce \
    fastq-dump "$source_url" --split-files -X "$spots" --outdir /data
  test -s "$out_dir/${run}_1.fastq"
  test -s "$out_dir/${run}_2.fastq"
  gzip -n -f "$out_dir/${run}_1.fastq" "$out_dir/${run}_2.fastq"
done
{
  printf 'artifact\tsha256\tsize_bytes\n'
  for f in "$out_dir"/*.fastq.gz; do printf '%s\t%s\t%s\n' "$(basename "$f")" "$(sha256sum "$f" | cut -d' ' -f1)" "$(stat -c%s "$f")"; done
} > "$out_dir/checksums.tsv"
