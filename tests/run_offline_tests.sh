#!/usr/bin/env bash
set -euo pipefail
python3 tests/create_synthetic_fixtures.py
python3 -m unittest discover -s tests -p 'test_*.py' -v
nextflow run main.nf -profile test --input tests/fixtures/generated/raw_samplesheet.tsv --fasta tests/fixtures/generated/reference.fa --gtf tests/fixtures/generated/reference.gtf --outdir tests/output/raw -resume
nextflow run main.nf -profile test --input tests/fixtures/generated/processed_samplesheet.tsv --fasta tests/fixtures/generated/reference.fa --gtf tests/fixtures/generated/reference.gtf --outdir tests/output/processed --run_de true --contrast condition,treatment,control -resume
nextflow run main.nf -profile test --input tests/fixtures/generated/one_vs_one.tsv --fasta tests/fixtures/generated/reference.fa --gtf tests/fixtures/generated/reference.gtf --outdir tests/output/guard_1v1 --run_de true --contrast condition,treatment,control -resume
nextflow run main.nf -profile test --input tests/fixtures/generated/processed_samplesheet.tsv --fasta tests/fixtures/generated/reference.fa --gtf tests/fixtures/generated/reference.gtf --outdir tests/output/guard_no_contrast --run_de true -resume
nextflow run main.nf -profile test --input tests/fixtures/generated/processed_samplesheet.tsv --fasta tests/fixtures/generated/reference.fa --gtf tests/fixtures/generated/reference.gtf --outdir tests/output/guard_rank --run_de true --design_formula '~ batch + condition' --contrast condition,treatment,control -resume
nextflow run main.nf -profile test --input tests/fixtures/generated/technical_replicates.tsv --fasta tests/fixtures/generated/reference.fa --gtf tests/fixtures/generated/reference.gtf --outdir tests/output/technical_replicates --run_de true --contrast condition,treatment,control -resume
python3 tests/assert_outputs.py
