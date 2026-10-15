#!/usr/bin/env bash
set -euo pipefail
out_dir=${1:-tests/integration/genuine_rna}
spots=${2:-100000}
source_url='https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR12413947/SRR12413947'
mkdir -p "$out_dir"
rm -f "$out_dir/SRR12413947_1.fastq.gz" "$out_dir/SRR12413947_2.fastq.gz" "$out_dir/SRR12413947.fastq.gz"
docker run --rm -v "$(pwd)/$out_dir:/data" -w /data \
  ncbi/sra-tools:3.4.1@sha256:4f048ba6754246e1ed7799a2c8ed469a38427f392e403ca08b93c4255104c9ce \
  fastq-dump "$source_url" --split-files -X "$spots" --outdir /data
if [[ -s "$out_dir/SRR12413947_1.fastq" && -s "$out_dir/SRR12413947_2.fastq" ]]; then
  gzip -n -f "$out_dir/SRR12413947_1.fastq" "$out_dir/SRR12413947_2.fastq"
elif [[ ! -e "$out_dir/SRR12413947_1.fastq" && -s "$out_dir/SRR12413947_2.fastq" ]]; then
  mv "$out_dir/SRR12413947_2.fastq" "$out_dir/SRR12413947.fastq"
  gzip -n -f "$out_dir/SRR12413947.fastq"
else
  printf 'Unexpected extracted read structure\n' >&2
  exit 1
fi
{
  printf 'artifact\tsha256\tsize_bytes\n'
  for f in "$out_dir"/*.fastq.gz; do
    printf '%s\t%s\t%s\n' "$(basename "$f")" "$(sha256sum "$f" | cut -d' ' -f1)" "$(stat -c%s "$f")"
  done
} > "$out_dir/checksums.tsv"
