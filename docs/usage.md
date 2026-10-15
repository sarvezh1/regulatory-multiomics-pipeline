# Usage

Run `nextflow run main.nf --pipeline_help true` for the compact help surface. Commands below use the validated parameter names. Replace logical paths with local inputs; explicit FASTA and GTF inputs are always required.

## Assay routes

RNA raw and processed inputs use the same assay selector; row-level raw or processed fields determine the route:

```bash
nextflow run main.nf -profile docker --assay RNA --input rna_samplesheet.tsv --fasta GRCh38.fa --gtf gencode.v50.gtf
nextflow run main.nf -profile docker --assay RNA --input rna_processed_samplesheet.tsv --fasta GRCh38.fa --gtf gencode.v50.gtf --run_de true --design_formula '~ condition' --contrast condition,treatment,control
```

ATAC raw and processed routes:

```bash
nextflow run main.nf -profile docker --assay ATAC --atac_input_route raw --input atac_samplesheet.tsv --fasta GRCh38.fa --gtf gencode.v50.gtf --blacklist GRCh38.blacklist.bed
nextflow run main.nf -profile docker --assay ATAC --atac_input_route processed --input atac_processed_samplesheet.tsv --fasta GRCh38.fa --gtf gencode.v50.gtf --blacklist GRCh38.blacklist.bed --run_da true --contrast condition,treatment,control
```

WGBS raw and processed routes:

```bash
nextflow run main.nf -profile docker --assay WGBS --wgbs_input_route raw --input wgbs_samplesheet.tsv --fasta GRCh38.fa --gtf gencode.v50.gtf
nextflow run main.nf -profile docker --assay WGBS --wgbs_input_route processed --input wgbs_processed_samplesheet.tsv --fasta GRCh38.fa --gtf gencode.v50.gtf --run_dm true --contrast condition,treatment,control
```

H3K27ac raw and processed routes:

```bash
nextflow run main.nf -profile docker --assay H3K27ac --chip_input_route raw --input chip_samplesheet.tsv --fasta GRCh38.fa --gtf gencode.v50.gtf --blacklist GRCh38.blacklist.bed
nextflow run main.nf -profile docker --assay H3K27ac --chip_input_route processed --input chip_processed_samplesheet.tsv --fasta GRCh38.fa --gtf gencode.v50.gtf --blacklist GRCh38.blacklist.bed --run_chip_diff true --contrast condition,treatment,control
```

Integration baseline, differential and partial-modality runs use the same manifest contract:

```bash
nextflow run main.nf -profile docker --assay INTEGRATION --integration_manifest integration_manifest.tsv --fasta GRCh38.fa --gtf gencode.v50.gtf --integration_mode baseline --integration_match_type REFERENCE_COMPARABLE
nextflow run main.nf -profile docker --assay INTEGRATION --integration_manifest integration_manifest.tsv --fasta GRCh38.fa --gtf gencode.v50.gtf --integration_mode differential --integration_match_type EXACT_SAMPLE_MATCHED
nextflow run main.nf -profile docker --assay INTEGRATION --integration_manifest partial_manifest.tsv --fasta GRCh38.fa --gtf gencode.v50.gtf --integration_mode both --integration_match_type CONDITION_MATCHED
```

Differential analyses require an explicit valid contrast, adequate independent biological replication and an estimable design. Technical libraries do not increase biological sample size.
