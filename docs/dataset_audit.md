# Dataset audit and validation ranking

Audit date: 2026-10-12. `VERIFIED` means the stated fact was confirmed in primary repository metadata or the primary publication. `PARTIALLY_VERIFIED` means core accessions/modalities were confirmed but matching, file-level detail, or unrestricted access remains incomplete. `UNVERIFIED` is not used as evidence. `REJECTED` means unsuitable as a primary v0.1 validation dataset, not scientifically invalid.

## ENCODE K562 reference epigenome

**ENCSR612NLL — VERIFIED as the strongest cross-assay reference/coordinate/contract validation resource.** The reference-epigenome record is released, describes human K562, lists 17 experiments, and includes all four anchors. The four experiments are separate records produced by different laboratories and biosample accessions. They are not the exact same biological material and must not be called sample-matched. The defensible description is **independent K562 experiments/cultures grouped into one reference epigenome**, not matched aliquots. Common cell-line identity and the ENCSR612NLL grouping support coordinate and contract validation, not differential multi-omics inference.

| Experiment | Assay | Biosamples / replicates | Platform/layout | Assembly and role | Status/audits |
|---|---|---|---|---|---|
| ENCSR109IQO | total RNA-seq | ENCBS538WJZ bio 1 tech 1; ENCBS864OKZ bio 2 tech 1 | HiSeq 2000; paired-end, reverse-stranded | released GRCh38 quantification/signal; also legacy hg19 | released; 1 warning |
| ENCSR868FGK | ATAC-seq | ENCBS049FLD bio 1 tech 1; ENCBS030HJI bio 2 tech 1; ENCBS254OQS bio 3 tech 1 | NovaSeq 6000; paired-end | released GRCh38 peaks/signal | released; 1 warning |
| ENCSR765JPC | WGBS | ENCBS316CKF bio 1 tech 1; ENCBS577JGM bio 2 tech 1 | HiSeq 2500; paired-end 100/125 nt | released GRCh38 methylation and alignments | released; 1 error, 1 non-compliance, 2 warnings |
| ENCSR000AKP | H3K27ac ChIP-seq | ENCBS639AAA bio 1 tech 1; ENCBS674MPN bio 2 tech 1; ENCBS378LYQ bio 3 tech 1 and 2 | single-end; platform record ENCSR000AKY | released GRCh38 and legacy hg19 outputs | released; 1 error, 9 warnings |

No donor relationship is applicable beyond K562 cell-line provenance; donor identifiers are not asserted. No treatment is stated for the four anchor experiments. H3K27ac control associations are present in the ENCODE record graph, but exact control accessions were not reliably exposed in the retrieved primary page and remain unresolved. The machine-readable audit lists verified representative raw/processed files and explicitly marks unresolved fields.

The current reference-epigenome view exposes compact selected GRCh38 outputs: RNA gene quantification/signal; ATAC peak/signal; WGBS CpG methylation tracks; and H3K27ac peaks/signal. These are suitable for Tier 3 processed integration. Raw WGBS alone includes individual FASTQs of roughly 40–74 GB and alignment files above 100 GB, demonstrating why it is unsuitable for local routine validation.

## GSE152146

**GSE152146 — VERIFIED at series and all 58 sample records; strongest biological study, with two separate validation roles.** The study uses HEK293 cells carrying a reporter construct and contains 38 RNA-seq, 8 ATAC-seq, 6 WGBS, 4 H3K27ac ChIP-seq, and 2 p65 ChIP-seq samples. Raw reads are in SRP266653/PRJNA638409 (58 SRA experiments; repository reports about 0.37 TB), and every sample has one processed file. The aggregate processed archive is 5.2 GB and was not downloaded.

1. **Four-modality baseline/integration case:** untreated RNA + ATAC + WGBS + H3K27ac. This role validates contract compatibility, coordinate integration, missingness handling, and baseline regulatory association. It is condition-comparable rather than proven replicate- or aliquot-matched.
2. **Partial-modality differential validation case:** only contrasts and modalities with legitimate independent biological replication may undergo differential inference. H3K27ac and WGBS observations without such replication remain descriptive/non-inferential by default.

GSE152146 is not a replicated four-modality differential benchmark.

The repository labels all reads paired-end: RNA on BGISEQ-500; ATAC, WGBS, H3K27ac, and p65 ChIP on NovaSeq 6000. Processing metadata declares hg38/GRCh38. Technical replicate numbers are not reported. Each GSM maps to a distinct BioSample and SRA experiment, so modalities are **comparable at condition level, not proven aliquots of the same biological replicate**.

### Target-modality condition matrix

Counts are biological libraries in repository metadata. Parentheses identify genotype where it differs from the reporter-parent/WT series.

| Condition | RNA | ATAC | WGBS | H3K27ac | Interpretation |
|---|---:|---:|---:|---:|---|
| untreated parent/control | 2 + a separate 2-sample dose cohort | 2 | 1 | 1 | only state containing all four; not replicate-matched across assays |
| 12 h TNF | 2 | 2 | 0 | 1 | RNA+ATAC replicated; H3K27ac descriptive |
| 12 d TNF | 2 | 0 | 1 | 0 | RNA+WGBS only; WGBS unreplicated |
| 12 d TNF + 10 d withdrawal | 2 + separate 2-sample dose cohort | 2 | 0 | 1 | RNA+ATAC replicated; H3K27ac descriptive |
| withdrawal + 12 h re-stimulation | 2 | 2 | 0 | 1 | RNA+ATAC replicated; H3K27ac descriptive |
| dose response, 0.4/50 ng/mL, naive and pretreated | 12 RNA libraries total including controls | 0 | 0 | 0 | RNA-only dose/memory design |
| WT untreated/9 h/recovery/re-stimulation | 8 | 0 | 0 | 0 | RNA-only genotype cohort |
| TET1/2/3 TKO same four states | 8 | 0 | 0 | 0 | RNA-only; differs from WGBS TET2-KO genotype |
| RELA KO untreated/12 d | 0 | 0 | 2 × 1 | 0 | WGBS-only, unreplicated per state |
| TET2 KO untreated/12 d | 0 | 0 | 2 × 1 | 0 | WGBS-only; not the RNA TET triple knockout |

There is **no untreated-versus-treatment contrast with all four modalities**. The strongest valid differential comparison in the target modalities is untreated versus 12 h TNF for RNA+ATAC (two biological libraries per group). Untreated versus recovery and recovery versus re-stimulation are also defensible RNA+ATAC contrasts. H3K27ac may be integrated descriptively at those conditions but cannot undergo valid replicate-based differential inference. Untreated versus 12 d has RNA+WGBS coverage, but WGBS is one library per condition and therefore descriptive under v0.1 rules.

The publication reports DMR results from these WGBS libraries; reproducing a published unreplicated comparison does not make it valid as a default inferential test. The full sample inventory is in `metadata/gse152146_sample_audit.tsv`.

## Additional candidates

### Human pancreatic islet cytokine response — PARTIALLY_VERIFIED

GSE108413 (RNA), GSE123404 (ATAC), and GSE133135 (H3K27ac), PMID 31676868, profile human islets treated with IFNγ + IL-1β for 48 h. A primary public data-package record reports RNA n=5, ATAC n=5, and H3K27ac n=4, with only partial donor overlap among assays. It has a strong inflammatory perturbation and three target modalities but no WGBS. Rank: best fallback biological dataset for RNA/ATAC/H3K27ac branch validation after sample-level matching is completed.

### NIH Roadmap Epigenomics — PARTIALLY_VERIFIED

The primary portal and consortium descriptions provide human reference epigenomes with RNA-seq, WGBS or related methylation assays, DNase accessibility, and H3K27ac. The processed portal reports 98 epigenomes with H3K27ac and 95 with a methylation assay; availability is broad and public. It generally uses DNase rather than ATAC, older hg19-era coordinates, and reference tissues are not differential matched cohorts. Rank: fallback reference/coordinate resource, not a primary v0.1 biological dataset.

### BLUEPRINT hematopoietic epigenomes — REJECTED as primary; processed fallback

BLUEPRINT generated RNA-seq, DNA methylation sequencing, DNase-seq, and H3K27ac among other ChIP assays in primary hematopoietic cells. Processed outputs are public, but the project states that most raw human data require managed EGA access. It lacks ATAC in the target combination and sample-level unrestricted matching was not completed here. Useful only as a processed primary-cell fallback.

### eGTEx — REJECTED as primary; future processed reference

eGTEx adds H3K27ac and WGBS to the extensive GTEx RNA resource and is scientifically valuable for donor-aware tissue regulation. GTEx states that raw DNA/RNA sequence and full donor metadata are protected access. ATAC is not part of the audited core combination. Public processed outputs may later support selected reference analyses, but access and scale make it unsuitable for v0.1 executable validation.

## Ranked strategy

1. **Technical/reference:** ENCSR612NLL, using selected released GRCh38 processed outputs. It has all four modalities and clear coordinate resources, but independent K562 cultures and no biological contrast.
2. **Biological baseline plus partial-modality differential:** GSE152146. Untreated provides the four-modality baseline/integration case; only predeclared replicated modality subsets provide differential validation. It is not a replicated four-modality differential benchmark.
3. **Fallback:** the pancreatic-islet cytokine datasets for replicated three-modality perturbation; Roadmap for reference methylation/H3K27ac/RNA; BLUEPRINT only for unrestricted processed files.

## Unresolved items

- Exact ENCODE audit-category text and all current downloadable raw-file rows could not be captured reliably from the dynamic portal; category counts and selected files are verified, and missing values are marked `UNRESOLVED`. Exact raw-file metadata must be resolved before raw H3K27ac validation.
- Exact per-file sizes for most GSE152146 supplementary files were not included in the family SOFT record. Filenames and availability are verified; unresolved sizes are not estimated.
- GSE152146 H3K27ac input-control identities are not represented as separate samples in the series metadata. They must be resolved before raw H3K27ac validation.
- Cross-assay aliquot lineage is not established for GSE152146; condition equivalence must not be upgraded to sample matching.

The two raw H3K27ac metadata gates above block raw H3K27ac validation only. They do not block processed H3K27ac descriptive integration or development and validation of RNA, ATAC, WGBS, shared contracts, or reporting.

## Primary sources

- [ENCODE ENCSR612NLL](https://www.encodeproject.org/reference-epigenomes/ENCSR612NLL/)
- [ENCODE ENCSR109IQO](https://www.encodeproject.org/experiments/ENCSR109IQO/), [ENCSR868FGK](https://www.encodeproject.org/experiments/ENCSR868FGK/), [ENCSR765JPC](https://www.encodeproject.org/experiments/ENCSR765JPC/), [ENCSR000AKP](https://www.encodeproject.org/experiments/ENCSR000AKP/)
- [GEO GSE152146](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152146), [BioProject PRJNA638409](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA638409), [primary publication, PMID 33155547](https://pmc.ncbi.nlm.nih.gov/articles/PMC7704108/)
- [Human-islet multi-omics package](https://data.pnnl.gov/group/nodes/dataset/34361)
- [NIH Roadmap GEO matrix](https://www.ncbi.nlm.nih.gov/geo/roadmap/epigenomics/?view=matrix), [processed-data portal](https://egg2.wustl.edu/roadmap/web_portal/processed_data)
- [BLUEPRINT data access](https://projects.ensembl.org/blueprint/)
- [GTEx data access](https://gtexportal.org/home/downloads/adult-gtex/)
