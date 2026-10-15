# Reproducibility validation

Reproducibility validation hardened reporting, provenance and replay behavior without changing assay or integration science. Validation used Nextflow 26.04.6, Docker Linux/amd64 and Python 3.12.7 for the host-side unified reporter. Existing pinned assay containers and tool versions remain authoritative in `assets/tool_containers.tsv` and generated `tool_manifest.tsv` / `container_manifest.tsv`.

## Unified reporting validation

Two fresh synthetic output/work trees were created independently for RNA, ATAC, WGBS, H3K27ac and full four-modality integration. Neither primary replay used a prior scientific-task cache. Each unified report contained 22 component statuses and 445 inventoried workflow artifacts before operational exclusions, plus schema-versioned run, input, reference, source, tool, container, parameter, output, status, warning, execution and plot manifests.

The normalized comparison found 17 selected canonical scientific/status/QC tables and classified them `IDENTICAL`. HTML, runtime records and operational reports were classified `SCIENTIFICALLY_EQUIVALENT_EXPECTED_RUNTIME_VARIATION`. Both reports produced the same stable run identity and source fingerprint when built from the same source, inputs, parameters and statuses.

Representative clean replay duration was approximately 9.5 minutes per five-route pass. Final report generation took approximately 1.2 seconds. Across 92 trace records in replay A, the observed maximum recorded task RSS was 1 GB. These are local synthetic execution observations, not production resource guarantees.

## Resume validation

An unchanged integration rerun with `-resume` reused all three scientific tasks. Changing the scientifically relevant distal maximum from 100,000 bp to 50,000 bp retained cached reference validation/build tasks and changed the regulatory-integration task hash from `bc/4acb7e` to `63/9d6b44`. Report generation occurs downstream of scientific execution and does not alter scientific task hashes.

## Failure and privacy guards

Automated reporting tests cover malformed TSV, malformed JSON, checksum mismatch, missing required scientific output behind a `COMPLETE` state, missing reference, private absolute path, duplicate output identity, source-fingerprint change and legitimate `NOT_RUN_DISABLED`. Representative report counts are checked directly against canonical input, output, status and warning manifests. The report generator rejects private-path or credential-like values and emits only relative logical paths.

Raw Nextflow trace, report, timeline and DAG files are retained as local operational metadata but excluded from the public output inventory. The sanitized execution table omits work paths and task directories.

## Upstream regression

- RNA: passed; 45 discovered contract/guard tests and all raw, processed, technical-replicate, differential and output assertions.
- ATAC: passed; 45 discovered contract/guard tests and all paired-end raw, processed differential, zero-peak and output assertions.
- WGBS: passed; 14 WGBS contract tests and nine raw/processed/guard workflow scenarios with output assertions.
- H3K27ac: passed; five ChIP contract tests and single-end, paired-end, processed differential and guard workflow assertions.
- Integration: passed; four unit tests and full, partial and contrast-mismatch workflow/output assertions.

The Nextflow update notice was observed but the validated 26.04.6 version was retained.

## Footprint and inclusion boundary

Retained reproducibility replay, report, resume and work evidence occupies 285,601,725 bytes (about 285.6 MB decimal; 272.4 MiB). Compact synthetic inputs are deterministic. Genuine FASTQ, BAM, bigWig, downloaded references, operational logs, work caches and generated results remain excluded from future repository inclusion. D: retained 89,414,295,552 bytes (about 89.4 GB decimal) free after replay and regression activity.

No new biological method, assay, regulatory link, threshold or interpretation was introduced.

## Pre-release status

Reproducibility validation is ready for release-stage review. An explicit project license and final citation metadata are not yet present and remain release-stage decisions. No repository, release, tag or DOI was created.
