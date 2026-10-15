# Reporting and run provenance

The unified reporting layer consumes completed workflow outputs. It does not rerun scientific analyses or create new biological inference.

`bin/build_run_report.py` creates a public-safe report directory containing the run, input, reference, source, tool, container, parameter, output, status, warning, execution and plot manifests defined in `schemas/run_provenance_contract.tsv`. Each tabular manifest carries `schema_name` and `schema_version`. `bin/validate_run_report.py` verifies required artifacts, manifest checksums, row structure, report-source counts, duplicate output identities and private-path hygiene.

When external-validation outputs are present, the report adds a separate external-validation section linking functional benchmark, annotation-concordance, biological-plausibility, status and provenance tables. The warning registry identifies cCRE dependence and the experimentally tested scope of the CRISPR benchmark. It never combines those validation roles into one score. External-validation SVG figures have same-stem TSV source tables and therefore enter the existing plot traceability manifest.

The run identifier is derived from normalized source, input, reference, tool, container, parameter and status identities. A timestamp is retained as operational context but is not the sole reproducibility identity. The source fingerprint hashes a sorted canonical list of execution-relevant files; generated outputs, work directories, caches, logs and downloaded data are excluded.

The offline HTML report separates process status, QC measurements and statistical inference. It links existing modality MultiQC reports rather than merging unlike QC systems. Integration results remain in the integration report and contracts. All displayed counts are sourced from `report_summary.tsv`, `status_manifest.tsv`, `output_manifest.tsv`, `warning_manifest.tsv` or `plot_manifest.tsv`.

Raw Nextflow traces, timelines, reports and DAG files are local operational records because they can contain local paths and task identifiers. They are excluded from the public output inventory. `sanitized_execution.tsv` retains selected process, status, duration, resource and container fields without task work paths.

Build and validate a report:

```bash
python bin/build_run_report.py \
  --results results \
  --source-root . \
  --output run_report \
  --requested-assay RNA,ATAC,WGBS,H3K27ac,INTEGRATION \
  --run-mode synthetic_validation \
  --nextflow-version 26.04.6 \
  --input-record samplesheet.tsv \
  --reference FASTA=GRCh38.primary_assembly.fa \
  --reference GTF=gencode.v50.annotation.gtf \
  --parameter-file effective_parameters.json \
  --tool-containers assets/tool_containers.tsv
python bin/validate_run_report.py run_report
```

The report generator is offline once workflow outputs, references and input metadata are present. It does not use remote scripts, fonts, analytics or report uploads.
