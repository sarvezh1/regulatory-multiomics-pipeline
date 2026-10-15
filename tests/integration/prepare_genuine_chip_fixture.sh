#!/usr/bin/env bash
set -euo pipefail
out_dir=${1:-tests/integration/genuine_chip}
reads=${2:-500000}
image='ncbi/sra-tools:3.4.1@sha256:4f048ba6754246e1ed7799a2c8ed469a38427f392e403ca08b93c4255104c9ce'
mkdir -p "$out_dir"
for run in SRR227385 SRR227650; do
  source_url="https://sra-pub-run-odp.s3.amazonaws.com/sra/${run}/${run}"
  rm -f "$out_dir/${run}.fastq" "$out_dir/${run}.fastq.gz"
  docker run --rm -v "$(pwd)/$out_dir:/data" -w /data "$image" \
    fastq-dump "$source_url" -X "$reads" --outdir /data
  test -s "$out_dir/${run}.fastq"
  gzip -n -f "$out_dir/${run}.fastq"
done
{
  printf 'artifact\tsha256\tsize_bytes\n'
  for file in "$out_dir"/*.fastq.gz; do
    printf '%s\t%s\t%s\n' "$(basename "$file")" "$(sha256sum "$file" | cut -d' ' -f1)" "$(stat -c%s "$file")"
  done
} > "$out_dir/checksums.tsv"
