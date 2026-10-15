#!/usr/bin/env python3
import argparse, csv, gzip, hashlib
from pathlib import Path

p=argparse.ArgumentParser(); p.add_argument("--input-dir",required=True); p.add_argument("--output-dir",required=True); a=p.parse_args()
src=Path(a.input_dir); out=Path(a.output_dir); out.mkdir(parents=True,exist_ok=True)
def write(path,fields,rows):
    with path.open("w",encoding="utf-8",newline="") as h:
        w=csv.DictWriter(h,fieldnames=fields,delimiter="\t",lineterminator="\n"); w.writeheader(); w.writerows(rows)
def digest(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1048576),b""): h.update(b)
    return h.hexdigest()

# ENCODE RSEM gene quantification includes non-Ensembl features; only original ENSG rows are safely harmonized.
rna=[]
with (src/"ENCFF611MXW.tsv").open(encoding="utf-8",newline="") as h:
    for r in csv.DictReader(h,delimiter="\t"):
        gid=r["gene_id"]
        if not gid.startswith("ENSG"): continue
        rna.append({"gene_id_original":gid,"gene_id_base":gid.rsplit(".",1)[0],"gene_name":"NA","normalized_expression":r["TPM"]})
write(out/"k562.rna_baseline.tsv",["gene_id_original","gene_id_base","gene_name","normalized_expression"],rna)

def peaks(accession,filename,prefix,outname):
    unique={}
    with gzip.open(src/filename,"rt",encoding="utf-8") as h:
        for line in h:
            c=line.rstrip("\n").split("\t")
            if c[0]!="chr22": continue
            start,end=int(c[1]),int(c[2]); token=hashlib.sha256(f"{c[0]}:{start}-{end}".encode()).hexdigest()[:16]
            unique[(start,end)]={"region_id":f"{prefix}_{token}","chromosome":c[0],"start":start,"end":end,"coordinate_system":"0-based-half-open"}
    rows=list(unique.values())
    rows.sort(key=lambda r:(r["start"],r["end"],r["region_id"]))
    write(out/outname,["region_id","chromosome","start","end","coordinate_system"],rows)
    return len(rows)
atac_n=peaks("ENCFF057UYP","ENCFF057UYP.bed.gz","ATAC_ENCFF057UYP","k562.atac_regions.tsv")
chip_n=peaks("ENCFF038DDS","ENCFF038DDS.bed.gz","H3K27AC_ENCFF038DDS","k562.h3k27ac_regions.tsv")

header=["modality","data_type","path","source_genome_build","coordinate_normalization_status","source_annotation","status","contrast_variable","numerator","denominator","integration_match_type"]
def row(mod,typ,name,annotation="NOT_APPLICABLE"):
    return {"modality":mod,"data_type":typ,"path":f"tests/integration/k562_integration/{name}","source_genome_build":"GRCh38","coordinate_normalization_status":"NORMALIZED" if mod!="RNA" else "NOT_APPLICABLE","source_annotation":annotation,"status":"COMPLETE","contrast_variable":"NA","numerator":"NA","denominator":"NA","integration_match_type":"REFERENCE_COMPARABLE"}
manifest=[row("RNA","rna_baseline","k562.rna_baseline.tsv","GENCODE V29 source quantification; identifiers retained"),row("ATAC","atac_regions","k562.atac_regions.tsv"),row("WGBS","wgbs_cpg_descriptive","k562.wgbs_chr22_descriptive.tsv"),row("H3K27ac","h3k27ac_regions","k562.h3k27ac_regions.tsv")]
write(out/"manifest.k562.tsv",header,manifest)
audit=[
 {"modality":"RNA","experiment":"ENCSR109IQO","file_accession":"ENCFF611MXW","file_type":"gene quantifications","assembly":"GRCh38","source_annotation":"V29","released":"true","original_size_bytes":10954745,"md5":"0b5bd054d75b1ec0df3b6a31f121d570","local_sha256":digest(src/"ENCFF611MXW.tsv"),"source_url":"https://www.encodeproject.org/files/ENCFF611MXW/","integration_use":"ENSG rows; TPM descriptive only"},
 {"modality":"ATAC","experiment":"ENCSR868FGK","file_accession":"ENCFF057UYP","file_type":"replicated narrowPeak","assembly":"GRCh38","source_annotation":"NOT_APPLICABLE","released":"true","original_size_bytes":6795339,"md5":"ae65def42f43bcf6d66b24b1509be56f","local_sha256":digest(src/"ENCFF057UYP.bed.gz"),"source_url":"https://www.encodeproject.org/files/ENCFF057UYP/","integration_use":"chr22 regions"},
 {"modality":"WGBS","experiment":"ENCSR765JPC","file_accession":"ENCFF750NZN","file_type":"bed9+ CpG bigBed","assembly":"GRCh38","source_annotation":"NOT_APPLICABLE","released":"true","original_size_bytes":1140822721,"md5":"9f055dab625d170a2392b7bffa00f227","local_sha256":"REMOTE_RANGE_READ","source_url":"https://www.encodeproject.org/files/ENCFF750NZN/","integration_use":"chr22 canonical CpGs; coverage plus rounded percent methylation; descriptive only"},
 {"modality":"H3K27ac","experiment":"ENCSR000AKP","file_accession":"ENCFF038DDS","file_type":"replicated narrowPeak","assembly":"GRCh38","source_annotation":"NOT_APPLICABLE","released":"true","original_size_bytes":1298162,"md5":"b24e283c2be89cbc2679eeb1688e4962","local_sha256":digest(src/"ENCFF038DDS.bed.gz"),"source_url":"https://www.encodeproject.org/files/ENCFF038DDS/","integration_use":"chr22 regions"},
]
write(out/"k562_source_audit.tsv",list(audit[0]),audit)
print(f"rna_genes={len(rna)} atac_chr22_peaks={atac_n} h3k27ac_chr22_peaks={chip_n}")
