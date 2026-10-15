# Reproducibility

## Validated boundary

The empirically validated local backend is Nextflow 26.04.6 with Docker on Linux/amd64 containers, invoked from Windows through WSL2. Other backends and architectures are not claimed as validated.

Scientific reproducibility depends on pinned input and reference identities, exact parameter snapshots, component statuses, executed tool versions, immutable container digests and the workflow-source fingerprint. The production reference policy remains GRCh38 primary assembly with GENCODE release 50 annotation context. Compact chromosome 22 and synthetic references are validation-only.

## Determinism

Canonical TSV contracts, genomic identifiers, coordinates, evidence classifications, component statuses and normalized provenance identities are expected to have deterministic scientific content for identical inputs and parameters. HTML files, logs, compressed files, binary formats and operational reports may contain timestamps, runtime fields, compression headers or tool-specific ordering. These are compared by normalized scientific content rather than unconditional byte identity.

The source fingerprint is the SHA-256 of a canonical sorted list of execution-relevant relative paths, file SHA-256 values and sizes. It includes `main.nf`, `nextflow.config`, `conf/`, `workflows/`, `modules/`, `bin/`, `schemas/` and relevant static assets. Work directories, results, genuine data, caches and generated provenance are excluded.

## Replay

Generate deterministic fixtures with the scripts under `tests/`, then run `tests/run_reproducibility_clean_replay.ps1`. The script uses fresh output and work directories for RNA, ATAC, WGBS, H3K27ac and integration and does not depend on an existing scientific-task cache. It performs a second fresh replay, compares normalized scientific outputs, and validates the unified reports.

Resume is tested separately with unchanged inputs and `-resume`. Nextflow decides cache reuse from process inputs, scripts, parameters and configuration. A scientifically relevant parameter change must invalidate affected tasks. Report generation runs after scientific tasks and does not invalidate them.

Bounded genuine validation procedures and accessions are documented in the block validation documents and `docs/public_data_sources.md`. Genuine sequencing data and downloaded reference indexes are not redistributed.

Exact byte reproduction can be limited by timestamps, runtime metadata, compressed headers and binary tool output. Such variation does not excuse changes in canonical scientific tables, deterministic identifiers, effect directions, missing-value semantics or component statuses.

## Hosted CI and local validation

GitHub CI performs public hygiene, documentation-link, release-metadata, TSV/JSON structure, contract, reporting, and external-validation unit checks. It has no genuine-data downloads and does not claim to reproduce container-heavy scientific validation.

The complete local validation record additionally includes Docker execution of all assay branches, genuine bounded public-data checks, two fresh synthetic replays, Nextflow cache-resume verification, and the full external K562 benchmark. This distinction keeps hosted checks practical without weakening the documented scientific evidence.
