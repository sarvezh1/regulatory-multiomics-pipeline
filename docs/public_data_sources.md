# Public validation data

The authoritative external-validation registry is `metadata/external_validation_sources.tsv`. It pins released K562 ATAC peaks `ENCFF057UYP`, released K562 H3K27ac peaks `ENCFF038DDS`, ENCODE Encyclopedia v4 K562 cCRE file `ENCFF286VQG` from `ENCSR913HQX`, and the K562 GRCh38 EPCrisprBenchmark training file from EngreitzLab/CRISPR_comparison commit `c075c1317a5d8bb47ad2fab33c8eb9d4bce99aac`. The registry includes retrieval date, checksum, role and independence class. Downloaded copies are validation-local and are not intended for repository inclusion.

The full-genome coordinate reference uses the separately pinned GENCODE release 50 GTF and GRCh38.p14 primary-assembly FASTA. These replace neither the production reference policy nor the compact chr22 execution fixtures.

Synthetic fixtures are generated deterministically by scripts under `tests/` and contain no private biological data. Bounded public raw subsets and downloaded references are local validation assets and are not intended for redistribution.

RNA validation used the bounded human run documented in `docs/validation/rna.md`. ATAC validation used the bounded paired-end human run documented in `docs/validation/atac.md`. WGBS validation used SRR4235744 / GSM2308597 as documented in `docs/validation/wgbs.md`. H3K27ac validation used ENCSR000AKP treatment ENCFF843VGI and matched control ENCFF245JEV as documented in `docs/validation/h3k27ac.md` and `metadata/encode_h3k27ac_control_mapping.tsv`.

Descriptive K562 integration uses the independent ENCODE anchors ENCSR109IQO (RNA), ENCSR868FGK (ATAC), ENCSR765JPC (WGBS) and ENCSR000AKP (H3K27ac). These experiments are `REFERENCE_COMPARABLE`, not matched aliquots, and support coordinate/contract validation rather than differential or causal claims. Exact local structural inputs are recorded in `metadata/k562_integration_sources.tsv` and `docs/validation/integration.md`.

Remote range reads, bounded extraction commands, source checksums and retained local checksums are documented in the relevant validation record. Full public experiments and publication-scale cohorts are not packaged with the project.
