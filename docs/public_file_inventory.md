# Public file inventory

This inventory defines what is intended for a future public repository. It does not itself initialize or stage a repository.

| Class | Inclusion policy | Examples |
|---|---|---|
| Public source | Include | `main.nf`, `nextflow.config`, `conf/`, `workflows/`, `modules/`, `bin/` |
| Public documentation | Include | `README.md`, `CHANGELOG.md`, `LICENSE`, `CITATION.cff`, `docs/` |
| Public schemas and metadata | Include | `schemas/`, curated tables under `metadata/`, `assets/tool_containers.tsv`, `assets/readme/` |
| Public small synthetic fixtures | Include selectively | deterministic fixture generators, compact fixture files needed by tests, and small retrieval/preparation scripts |
| Public validation summaries | Include | assay, integration, reproducibility, and external validation records; source registries; checksums; and compact summary tables where licensing permits |
| Generated reports | Exclude by default | run reports, workflow outputs, operational traces, timelines, and local DAGs |
| Local caches | Exclude | `work/`, `.nextflow*`, `tests/work/`, Python caches, container caches |
| Genuine or downloaded data | Exclude | FASTQ, BAM/BAI, bigWig/bigBed, SRA files, downloaded ENCODE and CRISPR resources |
| Downloaded references | Exclude | production FASTA/GTF, Bowtie2 indexes, Bismark references, blacklist downloads |
| Secrets and private state | Exclude | environment files, credentials, signed URLs, editor settings, local machine logs |

Project-authored code and documentation use the MIT License. Third-party data, annotations, containers, and software retain their original terms. Public source identifiers, retrieval instructions, and verified checksums are preferred over redistributing third-party datasets.

GitHub CI is intentionally lightweight: it checks public structure and hygiene, local documentation links, tabular/JSON syntax, release metadata, contract tests, reporting tests, and external-validation unit logic. Docker-heavy scientific workflows, clean-room replays, bounded genuine-data executions, and full external benchmark replays remain local validation and are documented rather than represented as hosted-CI jobs.
