#!/usr/bin/env python3
import csv
from pathlib import Path

ROOT=Path(__file__).resolve().parent/"fixtures"/"integration"
ROOT.mkdir(parents=True,exist_ok=True)
def tsv(name,fields,rows):
    p=ROOT/name; p.parent.mkdir(parents=True,exist_ok=True)
    with p.open("w",encoding="utf-8",newline="") as h:
        w=csv.DictWriter(h,fieldnames=fields,delimiter="\t",lineterminator="\n"); w.writeheader(); w.writerows(rows)
    return p

# Small GRCh38-compatible execution reference. Sequence content is immaterial to interval integration.
with (ROOT/"reference.fa").open("w",encoding="ascii",newline="\n") as h:
    h.write(">chr1\n")
    line="ACGT"*20
    for _ in range(12500): h.write(line+"\n")
genes=[
 ("ENSG00000000001.7","TX1","GENE1",3001,5000,"+"),("ENSG00000000002.2","TX2","GENE2",8001,10001,"-"),
 ("ENSG00000000003.1","TX3","GENE3",20001,22000,"+"),("ENSG00000000004.1","TX4","GENE4",30001,32000,"+"),
 ("ENSG00000000005.1","TX5","GENE5",40001,42000,"+"),("ENSG00000000006.1","TX6","GENE6",50001,52000,"+"),
 ("ENSG00000000007.1","TX7","GENE7",60001,62000,"+"),("ENSG00000000008.3","TX8A","GENE8",70001,72000,"+"),
 ("ENSG00000000008.3","TX8B","GENE8",75001,77000,"+"),("ENSG00000000009.1","TX9","BOUNDARY",2,1000,"+"),
 ("ENSG00000000010.1","TX10","TIE_A",200001,202000,"+"),("ENSG00000000011.1","TX11","TIE_B",208001,210001,"-"),
 ("ENSG00000000012.1","TX12","OVERLAP_A",90001,92000,"+"),("ENSG00000000013.1","TX13","OVERLAP_B",88051,90051,"-"),
 ("ENSG00000000001.7","TX1B","GENE1",3001,4800,"+"),
]
with (ROOT/"reference.gtf").open("w",encoding="utf-8",newline="\n") as h:
    h.write("##description: synthetic GENCODE v50-compatible integration fixture\n")
    for gid,tx,name,s,e,strand in genes:
        attrs=f'gene_id "{gid}"; transcript_id "{tx}"; gene_name "{name}";'
        h.write(f"chr1\tfixture\ttranscript\t{s}\t{e}\t.\t{strand}\t.\t{attrs}\n")

rna=[]
states={1:(2.5,1e-8),2:(-2.2,1e-7),3:(1.8,1e-5),4:(1.4,1e-4),5:(0.1,.8),6:(1.6,1e-5),7:(1.5,1e-5),8:(.2,.7),9:(0,.9),10:(0,.9),11:(0,.9),12:(0,.9),13:(0,.9)}
for i in range(1,14):
    gid=f"ENSG{i:011d}"; eff,p=states[i]
    rna.append({"gene_id_original":gid+((".7" if i==1 else ".3") if i in (1,8) else ".1"),"gene_id_base":gid,"gene_name":genes[i-1][2] if i<=len(genes) else f"GENE{i}","baseMean":100+i,"log2FoldChange":eff,"lfcSE":.2,"stat":eff/.2,"pvalue":p,"padj":p})
tsv("rna_differential.tsv",list(rna[0]),rna)
tsv("rna_baseline.tsv",["gene_id_original","gene_id_base","gene_name","normalized_expression"],[{"gene_id_original":r["gene_id_original"],"gene_id_base":r["gene_id_base"],"gene_name":r["gene_name"],"normalized_expression":100+i} for i,r in enumerate(rna)])

atac_specs=[
 ("ATAC_G1",2500,3300,2.0,1e-5),("ATAC_G1B",3100,3600,1.2,.01),("ATAC_G2",9500,10200,-2.0,1e-5),
 ("ATAC_G3",19500,20300,1.8,1e-4),("ATAC_G4",29500,30300,.1,.8),("ATAC_G5",39500,40300,2.0,1e-4),
 ("ATAC_G6",49500,50300,-1.7,1e-4),("ATAC_G7",59500,60300,1.5,1e-4),
 ("ATAC_G8A",69500,70300,1.3,.01),("ATAC_G8B",74500,75300,-1.4,.01),("ATAC_TWO_GENES",89900,90100,1.2,.01),
 ("ATAC_DISTAL",120000,120400,1.8,1e-4),("ATAC_TIE",205000,205001,1.2,.01),("ATAC_TOO_FAR",500000,500300,1.1,.01)
]
atac_regions=[{"region_id":x[0],"chromosome":"chr1","start":x[1],"end":x[2],"coordinate_system":"0-based-half-open"} for x in atac_specs]
atac_diff=[{"region_id":x[0],"chromosome":"chr1","start":x[1],"end":x[2],"baseMean":100,"log2FoldChange":x[3],"lfcSE":.2,"stat":x[3]/.2,"pvalue":x[4],"padj":x[4]} for x in atac_specs]
tsv("atac_regions.tsv",list(atac_regions[0]),atac_regions); tsv("atac_differential.tsv",list(atac_diff[0]),atac_diff)

chip_specs=[
 ("H3_G1",2600,3400,2.2,1e-5),("H3_G2",9600,10300,-2.1,1e-5),("H3_G4",29600,30400,.2,.8),
 ("H3_G5",39600,40400,2.0,1e-4),("H3_G7",59600,60400,1.7,1e-4),("H3_G8A",69600,70400,1.4,.01),
 ("H3_DISTAL",120100,120500,1.6,1e-4),("H3_TIE",204999,205002,1.1,.02)
]
chip_regions=[{"region_id":x[0],"chromosome":"chr1","start":x[1],"end":x[2],"coordinate_system":"0-based-half-open"} for x in chip_specs]
chip_diff=[{"region_id":x[0],"chromosome":"chr1","start":x[1],"end":x[2],"baseMean":120,"log2FoldChange":x[3],"lfcSE":.2,"stat":x[3]/.2,"pvalue":x[4],"padj":x[4]} for x in chip_specs]
tsv("h3k27ac_regions.tsv",list(chip_regions[0]),chip_regions); tsv("h3k27ac_differential.tsv",list(chip_diff[0]),chip_diff)

cpgs=[]
def add_cpg(gene_pos,m,u):
    for j,(sid,dm,du) in enumerate((("control_1",m,u),("control_2",m+1,u),("treatment_1",max(0,m-3),u+3),("treatment_2",max(0,m-2),u+2))):
        cov=dm+du; cpgs.append({"cpg_id":f"CpG_{gene_pos}_{j}","chromosome":"chr1","start":gene_pos,"end":gene_pos+2,"coordinate_system":"0-based-half-open","sample_id":sid,"methylated_count":dm,"unmethylated_count":du,"coverage":cov,"methylation_fraction":dm/cov,"measurement_state":"OBSERVED","passes_min_coverage":"true"})
for pos in (2900,3200,9800,20000,30000,40000,50000,60000,70000,75000,90000): add_cpg(pos,8,2)
tsv("wgbs_cpg.tsv",list(cpgs[0]),cpgs)
dmr_specs=[("DMR_G1",2800,3350,-.6,.01),("DMR_G2A",9550,9900,.5,.01),("DMR_G2B",9900,10300,.4,.02),("DMR_G4",29800,30200,0,.9),("DMR_G7",59800,60200,.5,.01),("DMR_G8A",69800,70200,-.4,.02),("DMR_G8B",74800,75200,.4,.02)]
dmrs=[{"dmr_id":x[0],"chromosome":"chr1","start":x[1],"end":x[2],"coordinate_system":"0-based-half-open","n_cpg":5,"mean_methylation_numerator":.3,"mean_methylation_denominator":.7,"delta_methylation":x[3],"delta_methylation_pp":x[3]*100,"test_statistic":4,"pvalue":x[4],"qvalue":x[4],"contrast":"condition,treatment,control","method":"dmrseq"} for x in dmr_specs]
tsv("wgbs_dmr.tsv",list(dmrs[0]),dmrs)

external=[{"source":"fixture_curated","version":"1","genome_build":"GRCh38","chromosome":"chr1","start":120000,"end":120400,"gene_identifier":"ENSG00000000005.1","evidence_type":"CRISPR","validation_type":"EXPERIMENTALLY_VALIDATED"}]
tsv("external_links.tsv",list(external[0]),external)

header=["modality","data_type","path","source_genome_build","coordinate_normalization_status","source_annotation","status","contrast_variable","numerator","denominator","integration_match_type"]
def row(mod,typ,name,status="COMPLETE",num="treatment",den="control"):
    return {"modality":mod,"data_type":typ,"path":f"tests/fixtures/integration/{name}","source_genome_build":"synthetic_GRCh38_compatible","coordinate_normalization_status":"NORMALIZED" if mod!="RNA" else "NOT_APPLICABLE","source_annotation":"GENCODE_50_fixture" if mod=="RNA" else "NOT_APPLICABLE","status":status,"contrast_variable":"condition" if "differential" in typ or typ=="wgbs_dmr" else "NA","numerator":num if "differential" in typ or typ=="wgbs_dmr" else "NA","denominator":den if "differential" in typ or typ=="wgbs_dmr" else "NA","integration_match_type":"EXACT_SAMPLE_MATCHED"}
full=[row("RNA","rna_baseline","rna_baseline.tsv"),row("RNA","rna_differential","rna_differential.tsv"),row("ATAC","atac_regions","atac_regions.tsv"),row("ATAC","atac_differential","atac_differential.tsv"),row("WGBS","wgbs_cpg","wgbs_cpg.tsv"),row("WGBS","wgbs_dmr","wgbs_dmr.tsv"),row("H3K27ac","h3k27ac_regions","h3k27ac_regions.tsv"),row("H3K27ac","h3k27ac_differential","h3k27ac_differential.tsv"),row("EXTERNAL","enhancer_gene_map","external_links.tsv",num="NA",den="NA")]
full[-1].update({"contrast_variable":"NA","numerator":"NA","denominator":"NA","source_genome_build":"GRCh38","coordinate_normalization_status":"NORMALIZED","source_annotation":"NOT_APPLICABLE"})
tsv("manifest.full.tsv",header,full)
tsv("manifest.partial.tsv",header,[r for r in full if r["modality"] in {"RNA","ATAC"}])
mismatch=[dict(r) for r in full]
for r in mismatch:
    if r["data_type"]=="atac_differential": r["numerator"]="recovery"
tsv("manifest.mismatch.tsv",header,mismatch)
print(ROOT)
