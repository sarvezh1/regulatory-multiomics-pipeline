#!/usr/bin/env python3
import argparse, csv, json
from pathlib import Path

def rows(path):
    with open(path,encoding="utf-8",newline="") as h: return list(csv.DictReader(h,delimiter="\t"))
def one(items,**want):
    hit=[r for r in items if all(r.get(k)==str(v) for k,v in want.items())]
    assert len(hit)==1,(want,len(hit)); return hit[0]

p=argparse.ArgumentParser(); p.add_argument("--full",required=True); p.add_argument("--partial",required=True); p.add_argument("--mismatch",required=True); a=p.parse_args()
full=Path(a.full)/"integration"; partial=Path(a.partial)/"integration"; mismatch=Path(a.mismatch)/"integration"
prom=rows(full/"reference/integration.promoters.m2000_p500.tsv")
tss=rows(full/"reference/integration.tss.tsv")
ev=rows(full/"evidence/promoter_evidence.tsv")
assert one(prom,gene_name="BOUNDARY")["promoter_start"]=="0"
assert len([r for r in prom if r["gene_name"]=="GENE8"])==2
assert "TX1;TX1B"==one(tss,gene_name="GENE1")["transcript_ids"]
g1=one(ev,gene_name="GENE1")
assert g1["gene_id_original"]=="ENSG00000000001.7"
assert g1["gene_id_base"]=="ENSG00000000001"
assert g1["regulatory_pattern"]=="RNA_UP__ATAC_UP__H3K27AC_UP__METHYLATION_HYPOMETHYLATED"
g2=one(ev,gene_name="GENE2")
assert g2["regulatory_pattern"]=="RNA_DOWN__ATAC_DOWN__H3K27AC_DOWN__METHYLATION_HYPERMETHYLATED"
g4=one(ev,gene_name="GENE4")
assert g4["rna_state"]=="UP"
assert g4["atac_state"]=="NO_SIGNIFICANT_CHANGE" and g4["h3k27ac_state"]=="NO_SIGNIFICANT_CHANGE"
g5=one(ev,gene_name="GENE5")
assert g5["rna_state"]=="NO_SIGNIFICANT_CHANGE"
assert g5["atac_state"]=="UP" and g5["h3k27ac_state"]=="UP"
assert int(g1["n_atac_regions"])==2
assert int(g2["n_dmr"])==2
assert one(ev,gene_name="GENE6")["regulatory_pattern"].startswith("RNA_UP__ATAC_DOWN")
assert one(ev,gene_name="GENE7")["regulatory_pattern"].endswith("METHYLATION_HYPERMETHYLATED")
gene8=one(rows(full/"evidence/gene_summary.tsv"),gene_name="GENE8")
assert gene8["number_of_promoters"]=="2" and gene8["promoter_interpretation"]=="PROMOTER_HETEROGENEITY"
two=[r for r in rows(full/"relationships/atac_promoter.tsv") if r["source_region_id"]=="ATAC_TWO_GENES"]
assert len({r["gene_id_base"] for r in two})==2
ties=[r for r in rows(full/"relationships/distal_candidate_links.tsv") if r["region_id"]=="ATAC_TIE"]
assert len(ties)==2 and {r["candidate_gene_id"] for r in ties}=={"ENSG00000000010","ENSG00000000011"}
assert not any(r["region_id"]=="ATAC_TOO_FAR" for r in rows(full/"relationships/distal_candidate_links.tsv"))
ext=rows(full/"relationships/external_enhancer_gene_links.tsv")
assert ext[0]["link_method"]=="EXTERNAL_EXPERIMENTALLY_VALIDATED"
coloc=rows(full/"relationships/atac_h3k27ac_colocalization.tsv")
assert any(r["atac_region_id"]=="ATAC_DISTAL" and r["relationship"]=="ATAC_H3K27AC_COLOCATED" for r in coloc)
meth=rows(full/"evidence/promoter_methylation_by_sample.tsv")
for r in meth:
    m=int(r["total_methylated_reads"]); u=int(r["total_unmethylated_reads"])
    assert abs(float(r["coverage_weighted_methylation_fraction"])-m/(m+u))<1e-12
sens=rows(full/"sensitivity/integration.sensitivity_summary.tsv")
assert {r["setting"] for r in sens if r["analysis"]=="PROMOTER_WINDOW"}=={"-2000/+500","m1000_p500","m2000_p1000"}
assert all(r["pattern_assignments_changed"]!="NA" for r in sens if r["analysis"]=="PROMOTER_WINDOW")
assert {r["setting"] for r in sens if r["analysis"]=="DISTAL_MAX_DISTANCE"}=={"50000","100000","250000"}
status={r["component"]:r["status"] for r in rows(full/"provenance/integration.analysis_status.tsv")}
assert status["integration_differential"]=="COMPLETE"
assert json.loads((full/"provenance/integration.provenance.json").read_text())["no_cross_assay_pvalue_combination"] is True
pev=rows(partial/"evidence/promoter_evidence.tsv")
assert all(r["h3k27ac_state"]=="NOT_AVAILABLE" and r["methylation_state"]=="NOT_AVAILABLE" for r in pev)
partial_g3=one(pev,gene_name="GENE3")
assert partial_g3["regulatory_pattern"]=="RNA_UP__ATAC_UP__H3K27AC_NOT_AVAILABLE__METHYLATION_NOT_AVAILABLE"
assert one(ev,gene_name="GENE3")["h3k27ac_state"]=="NOT_TESTED"
mstat={r["component"]:r["status"] for r in rows(mismatch/"provenance/integration.analysis_status.tsv")}
assert mstat["integration_differential"]=="NOT_RUN_INVALID_DESIGN"
assert all(r["compatibility"]=="CONTRAST_MISMATCH" for r in rows(mismatch/"compatibility/integration.contrast_compatibility.tsv"))
assert all(r["integration_match_type"]=="EXACT_SAMPLE_MATCHED" for r in ev)
print("Integration output assertions passed")
