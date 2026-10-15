# External validation record

External validation completed functional benchmarking, non-independent annotation concordance, and reference-comparable biological plausibility analysis without changing previously validated scientific outputs.

## Frozen inputs and full-genome predictions

The full validation reference contains 389,243 distinct GENCODE v50 TSS records on primary chromosomes. Released K562 inputs contributed 151,399 ATAC regions (`ENCFF057UYP`) and 52,318 H3K27ac regions (`ENCFF038DDS`). The validation construction produced 285,903 promoter–ATAC links, 264,274 promoter–H3K27ac links, 46,940 ATAC/H3K27ac overlaps, and 99,886 / 103,288 / 103,919 candidate-link rows at 50 / 100 / 250 kb.

The frozen 100 kb candidate-link SHA-256 is `98cd221f10559e95ce169a9866f05a84e853e3ac95a9ca29e4de959b54b79fb7`. The 50 kb and 250 kb checksums are recorded beside it in `metadata/external_validation_prediction_freeze.tsv`. The truth resource had not been loaded when these values were frozen.

## cCRE concordance

The current ENCODE Encyclopedia v4 K562 cCRE source is experiment `ENCSR913HQX`, file `ENCFF286VQG`, MD5 `344a9e65f749a3e5703de96a3577101b`, SHA-256 `f72211065ed101e9998a39676295c1e96e1abd2c21aed2025dc6e3a8100c01e0`. The file contains 2,348,854 primary cCRE records.

Any-cCRE overlap was 16,223/16,452 (98.61%) for distal colocated regions, 81,346/87,387 (93.09%) for distal ATAC-only regions, 5,944/7,187 (82.70%) for distal H3K27ac-only regions, 46,184/47,560 (97.11%) for promoter ATAC regions, and 27,797/28,327 (98.13%) for promoter H3K27ac regions. Exact source classes and overlaps are in `ccre_concordance_summary.tsv` and `ccre_region_overlaps.tsv`. These are annotation-concordance results, not an independent benchmark.

## Functional benchmark

The pinned CRISPR resource SHA-256 is `9eddfe18f6797ba33db8ba7e75378cc6af6b79cbb83f62d470c6e0cefed79405`. It contains 10,356 K562 GRCh38 pairs: 471 regulated positives and 9,885 negatives before project-specific eligibility filters. Source studies are Gasperini2019, Nasser2021, and Schraivogel2020.

Five rows had unresolved GENCODE v50 identifiers and 63 were promoter-proximal tests. The primary distal scoring denominator was therefore 10,288 pairs: 427 positives and 9,861 negatives.

At the frozen 100 kb setting:

| Metric | Value |
|---|---:|
| TP / FP / TN / FN | 71 / 51 / 9,810 / 356 |
| Precision | 0.5820 |
| Recall | 0.1663 |
| Specificity | 0.9948 |
| F1 | 0.2587 |
| Balanced accuracy | 0.5806 |
| MCC | 0.2968 |
| Fisher odds ratio | 38.36 |
| Fisher two-sided p-value | 2.56e-67 |
| Approximate odds-ratio 95% interval | 26.37–55.82 |

The predefined distance sensitivities produced 68/51/9,810/359 at 50 kb and 71/51/9,810/356 at 250 kb (TP/FP/TN/FN). The v0.1 production default remains 100 kb; these results did not retune it.

The 100 kb distance-only baseline produced TP 87, FP 62, TN 9,799 and FN 340 (F1 0.3021, MCC 0.3297). Thus, in this benchmark, the frozen activity-region requirement did not outperform simple distance-only nearest-TSS assignment on F1 or MCC. Activity-region overlap without gene assignment recovered 421/427 positives but also overlapped 8,994/9,861 negatives, showing why region detection and gene linking must be evaluated separately.

Restricting predicted positives to ATAC/H3K27ac-colocated candidate regions produced TP 67, FP 34, TN 9,827 and FN 360: precision 0.6634, recall 0.1569, F1 0.2538 and MCC 0.3104. This is descriptive and does not create a new production filter.

At 100 kb, source-study F1 values were 0.2910 for Gasperini2019 (n=5,248), 0.1538 for Nasser2021 (n=3,735), and 0.0800 for Schraivogel2020 (n=1,305). Experimental definitions differ across source studies, so aggregate results do not erase study provenance.

The primary failure breakdown was: 71 recovered positives; 51 candidate overlaps of experimentally negative pairs; 9,810 true negatives; and 356 false negatives. The false negatives comprise 6 elements without an activity region, 259 cases where the tested gene was not nearest, 78 beyond 100 kb, and 13 residual link non-recoveries. Row-level tables make every classification inspectable.

## Biological plausibility

The chr22 four-modality context includes 1,324 genes. Promoter ATAC evidence increased from 39.89% in the bottom expression stratum to 78.67% in the top stratum; promoter H3K27ac evidence increased from 9.04% to 66.83%. Median coverage-weighted promoter methylation was 0.1292, 0.2450 and 0.2611 across bottom, middle and top strata. These independent experiments are not matched aliquots, and the methylation tendency is not interpreted causally.

Descriptive cross-gene Spearman rho was 0.2629 for RNA versus promoter ATAC (n=1,324), 0.6102 for RNA versus promoter H3K27ac (n=1,324), and 0.1971 for RNA versus promoter methylation (n=1,310). P-values were intentionally not calculated.

## Determinism and status

The functional benchmark was run twice from the same frozen inputs. All seven canonical pair, exclusion, denominator, metric, distance-stratum and failure tables were byte-identical. Candidate-link checksums remained equal to the pre-comparison freeze. `external_source_audit`, `ccre_concordance`, `crispr_functional_benchmark`, and `biological_plausibility` are `COMPLETE`; optional ABC is `NOT_RUN_DISABLED`; optional compact genuine differential validation is `NOT_AVAILABLE_INPUT` for the recorded reason.

Python 3.12.7 and standard-library deterministic interval/statistical utilities were executed. No new production tool, assay, model, aligner index, raw cohort, ABC computation, Hi-C analysis, machine learning, or production parameter was introduced.

Measured on the validated local Windows/Docker host, full-genome candidate construction took 17.0 seconds with an observed peak working set of 1.54 GB; one CRISPR scoring replay took 48.9 seconds with 1.13 GB; and full cCRE concordance took 10.5 seconds with 339 MB. The retained external-validation reference, source, derived and report evidence occupies 1,516,232,468 bytes (about 1.52 GB decimal). D: retained 87,733,682,176 bytes (about 87.73 GB decimal) free at closure.

Two implementation issues were corrected without changing scientific rules. A timed-out range-transfer process had left an oversized partial GTF; it was stopped before scientific use and the official file was fetched cleanly and checksum-verified. The containerized integration reference builder was impractical for the full GENCODE annotation, so a streaming builder was introduced and validated field-for-field against every chr22 integration TSS/promoter record before full-genome use.

## Claim boundary

The benchmark supports a limited external evaluation of candidate-link recovery among experimentally tested K562 pairs. It does not establish genome-wide truth, causal methylation, matched multi-omics effects, physical contact, or that a nearest TSS is the true target. cCRE overlap remains non-independent annotation concordance. Docker Linux/amd64 remains the empirically validated production backend.
