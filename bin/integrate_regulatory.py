#!/usr/bin/env python3
import argparse, bisect, csv, hashlib, json, math, os, re
from collections import defaultdict
from pathlib import Path

NA="NA"; ALLOWED={f"chr{i}" for i in range(1,23)}|{"chrX","chrY"}
REQ_MANIFEST={"modality","data_type","path","source_genome_build","coordinate_normalization_status","source_annotation","status","contrast_variable","numerator","denominator","integration_match_type"}
TYPES={"RNA":{"rna_baseline","rna_differential"},"ATAC":{"atac_regions","atac_differential"},"WGBS":{"wgbs_cpg","wgbs_cpg_descriptive","wgbs_dmr"},"H3K27ac":{"h3k27ac_regions","h3k27ac_differential"},"EXTERNAL":{"enhancer_gene_map"}}

def read_tsv(path):
    with open(path,encoding="utf-8",newline="") as h: return list(csv.DictReader(h,delimiter="\t"))
def write_tsv(path, rows, fields):
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    with open(path,"w",encoding="utf-8",newline="") as h:
        w=csv.DictWriter(h,delimiter="\t",fieldnames=fields,lineterminator="\n",extrasaction="ignore"); w.writeheader()
        for r in rows: w.writerow({k:(NA if r.get(k) is None else r.get(k,NA)) for k in fields})
def sha(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for b in iter(lambda:f.read(1048576),b""): h.update(b)
    return h.hexdigest()
def fnum(x):
    if x in (None,"",NA,"NaN","nan"): return None
    try: return float(x)
    except ValueError: return None
def fint(x):
    v=fnum(x); return None if v is None else int(v)
def overlap(a0,a1,b0,b1): return max(0,min(a1,b1)-max(a0,b0))
def base_gene(x): return re.sub(r"\.[0-9]+$","",x or "")
def state(effect,padj,complete,fdr,methyl=False):
    if not complete: return "NOT_TESTED"
    e=fnum(effect); q=fnum(padj)
    if e is None or q is None: return "NOT_TESTED"
    if q>fdr or e==0: return "NO_SIGNIFICANT_CHANGE"
    if methyl: return "HYPERMETHYLATED" if e>0 else "HYPOMETHYLATED"
    return "UP" if e>0 else "DOWN"
def safe_join(values): return ";".join(sorted(set(str(x) for x in values if x not in (None,"",NA)))) or NA
def validate_interval_rows(rows,kind):
    seen=set()
    for i,r in enumerate(rows,2):
        chrom=r.get("chromosome"); s=fint(r.get("start")); e=fint(r.get("end"))
        if chrom not in ALLOWED or s is None or e is None or s<0 or e<=s: raise ValueError(f"{kind} invalid interval at row {i}")
        rid=r.get("region_id") or r.get("dmr_id") or r.get("cpg_id")
        if rid and (rid,chrom,s,e) in seen: raise ValueError(f"{kind} duplicate identifier/interval at row {i}")
        seen.add((rid,chrom,s,e))

def interval_index(rows,start_field="start",end_field="end"):
    grouped=defaultdict(list); maximum=0
    for r in rows:
        s=int(r[start_field]); e=int(r[end_field]); grouped[r["chromosome"]].append((s,e,r)); maximum=max(maximum,e-s)
    result={}
    for chrom,items in grouped.items():
        items.sort(key=lambda x:(x[0],x[1])); result[chrom]=([x[0] for x in items],items)
    return result,maximum

def indexed_hits(index,maximum_length,chrom,start,end):
    if chrom not in index: return []
    starts,items=index[chrom]
    lo=bisect.bisect_left(starts,start-maximum_length)
    hi=bisect.bisect_left(starts,end)
    return [r for s,e,r in items[lo:hi] if e>start]

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--manifest",required=True); ap.add_argument("--promoters",required=True); ap.add_argument("--tss",required=True); ap.add_argument("--sensitivity-promoters",nargs="*",default=[]); ap.add_argument("--outdir",required=True); ap.add_argument("--mode",choices=["baseline","differential","both"],default="both"); ap.add_argument("--fdr",type=float,default=.05); ap.add_argument("--distal-max",type=int,default=100000); ap.add_argument("--external-map"); ap.add_argument("--match-type",required=True,choices=["EXACT_SAMPLE_MATCHED","CONDITION_MATCHED","REFERENCE_COMPARABLE","UNMATCHED"]); a=ap.parse_args()
    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    manifest=read_tsv(a.manifest)
    if not manifest or not REQ_MANIFEST.issubset(manifest[0]): raise ValueError("integration manifest missing required columns")
    if {r["integration_match_type"] for r in manifest} != {a.match_type}:
        raise ValueError("manifest integration_match_type does not match the requested run class")
    loaded={}; compatibility=[]; checksums=[]
    for r in manifest:
        mod=r["modality"]; typ=r["data_type"]
        if mod not in TYPES or typ not in TYPES[mod]: raise ValueError(f"unsupported manifest modality/data_type: {mod}/{typ}")
        p=Path(r["path"])
        if not p.is_file(): raise ValueError(f"missing integration input: {p}")
        ref_ok=r["source_genome_build"] in {"GRCh38","GRCh38_primary","GRCh38.p14_primary","synthetic_GRCh38_compatible"}
        norm_ok=r["coordinate_normalization_status"] in {"COMPLETE","NORMALIZED","verified_primary_assembly","NOT_APPLICABLE"}
        coord_type=typ not in {"rna_baseline","rna_differential"}
        reference_state="REFERENCE_COMPATIBLE" if ref_ok and (norm_ok or not coord_type) else "REFERENCE_MISMATCH"
        compatibility.append({"modality":mod,"data_type":typ,"reference_compatibility":reference_state,"status":r["status"],"source_genome_build":r["source_genome_build"],"coordinate_normalization_status":r["coordinate_normalization_status"]})
        if reference_state!="REFERENCE_COMPATIBLE": raise ValueError(f"reference/coordinate compatibility failed for {mod}/{typ}")
        loaded[(mod,typ)]=read_tsv(p); checksums.append({"modality":mod,"data_type":typ,"filename":p.name,"sha256":sha(p)})
    write_tsv(out/"compatibility"/"integration.contract_compatibility.tsv",compatibility,["modality","data_type","reference_compatibility","status","source_genome_build","coordinate_normalization_status"])
    promoters=read_tsv(a.promoters); tss=read_tsv(a.tss)
    for p in promoters:
        p["start"]=int(p["promoter_start"]); p["end"]=int(p["promoter_end"])
    promoter_index,promoter_max_length=interval_index(promoters)
    for typ in ("atac_regions","atac_differential","h3k27ac_regions","h3k27ac_differential","wgbs_cpg","wgbs_cpg_descriptive","wgbs_dmr"):
        for mod in TYPES:
            if (mod,typ) in loaded: validate_interval_rows(loaded[(mod,typ)],typ)
    available={m:any(k[0]==m for k in loaded) for m in ("RNA","ATAC","WGBS","H3K27ac")}
    # Differential contrast compatibility is manifest-driven and never inferred from labels.
    diff_types={"rna_differential","atac_differential","wgbs_dmr","h3k27ac_differential"}
    diff_rows=[r for r in manifest if r["data_type"] in diff_types]
    triples={(r["contrast_variable"],r["numerator"],r["denominator"]) for r in diff_rows if r["status"]=="COMPLETE"}
    if not diff_rows or not triples: cstate="CONTRAST_UNAVAILABLE"
    elif len(triples)==1: cstate="CONTRAST_COMPATIBLE"
    else: cstate="CONTRAST_MISMATCH"
    contrast_rows=[{"modality":r["modality"],"data_type":r["data_type"],"contrast_variable":r["contrast_variable"],"numerator":r["numerator"],"denominator":r["denominator"],"upstream_status":r["status"],"compatibility":cstate} for r in diff_rows]
    write_tsv(out/"compatibility"/"integration.contrast_compatibility.tsv",contrast_rows,["modality","data_type","contrast_variable","numerator","denominator","upstream_status","compatibility"])
    differential_allowed=cstate=="CONTRAST_COMPATIBLE" and a.mode in {"differential","both"}
    # Gene-level RNA evidence.
    rna_base={base_gene(r.get("gene_id_base") or r.get("gene_id_original")):r for r in loaded.get(("RNA","rna_baseline"),[])}
    rna_diff={base_gene(r.get("gene_id_base") or r.get("gene_id_original")):r for r in loaded.get(("RNA","rna_differential"),[])}
    mapping=[]
    for typ,rows in (("rna_baseline",loaded.get(("RNA","rna_baseline"),[])),("rna_differential",loaded.get(("RNA","rna_differential"),[]))):
        for r in rows:
            source=r.get("gene_id_original",""); mapped=base_gene(r.get("gene_id_base") or source)
            mapping.append({"data_type":typ,"source_identifier":source,"mapped_identifier":mapped,"mapping_method":"REMOVE_TERMINAL_NUMERIC_ENSEMBL_VERSION","mapping_status":"MAPPED" if mapped.startswith("ENSG") else "UNRESOLVED"})
    write_tsv(out/"relationships"/"gene_identifier_mapping.tsv",mapping,["data_type","source_identifier","mapped_identifier","mapping_method","mapping_status"])
    diff_complete={r["modality"]:(r["status"]=="COMPLETE" and differential_allowed) for r in diff_rows}
    # Region/promoter detailed relationships.
    def link_regions(rows,idfield,prefix):
        links=[]
        for r in rows:
            rs,re_=int(r["start"]),int(r["end"])
            for p in indexed_hits(promoter_index,promoter_max_length,r["chromosome"],rs,re_):
                ov=overlap(rs,re_,p["start"],p["end"])
                if ov:
                    links.append({"source_region_id":r[idfield],"promoter_id":p["promoter_id"],"gene_id_base":p["gene_id_base"],"chromosome":r["chromosome"],"overlap_bp":ov,"source_region_length":re_-rs,"promoter_length":p["end"]-p["start"]})
        return links
    atac_regions=loaded.get(("ATAC","atac_regions"),[]); chip_regions=loaded.get(("H3K27ac","h3k27ac_regions"),[])
    atac_index,atac_max_length=interval_index(atac_regions) if atac_regions else ({},0)
    atac_links=link_regions(atac_regions,"region_id","ATAC"); chip_links=link_regions(chip_regions,"region_id","H3K27ac")
    write_tsv(out/"relationships"/"atac_promoter.tsv",atac_links,["source_region_id","promoter_id","gene_id_base","chromosome","overlap_bp","source_region_length","promoter_length"])
    write_tsv(out/"relationships"/"h3k27ac_promoter.tsv",chip_links,["source_region_id","promoter_id","gene_id_base","chromosome","overlap_bp","source_region_length","promoter_length"])
    # CpG links and count-weighted promoter/sample methylation.
    cpgs=loaded.get(("WGBS","wgbs_cpg"),[]); descriptive_cpgs=loaded.get(("WGBS","wgbs_cpg_descriptive"),[]); cpg_links=[]; meth=defaultdict(lambda:[0,0,set(),0]); meth_desc=defaultdict(lambda:[0.0,0,set()])
    for r in cpgs:
        s,e=int(r["start"]),int(r["end"])
        m=fint(r.get("methylated_count")); u=fint(r.get("unmethylated_count")); cov=fint(r.get("coverage"))
        if m is None or u is None or m<0 or u<0 or cov!=m+u: raise ValueError("invalid WGBS count identity")
        for p in indexed_hits(promoter_index,promoter_max_length,r["chromosome"],s,e):
            if overlap(s,e,p["start"],p["end"]):
                ov=overlap(s,e,p["start"],p["end"]); sid=r.get("sample_id",NA)
                cpg_links.append({"cpg_id":r["cpg_id"],"promoter_id":p["promoter_id"],"sample_id":sid,"overlap_bp":ov,"methylated_count":m,"unmethylated_count":u,"coverage":cov,"methylation_fraction":(m/cov if cov else NA),"evidence_type":"COUNT_PRESERVING"})
                z=meth[(p["promoter_id"],sid)]; z[0]+=m; z[1]+=u; z[2].add(r["cpg_id"]); z[3]+=1 if str(r.get("passes_min_coverage","")).lower()=="true" else 0
    for r in descriptive_cpgs:
        s,e=int(r["start"]),int(r["end"]); cov=fint(r.get("coverage")); frac=fnum(r.get("methylation_fraction"))
        if cov is None or cov<0 or frac is None or not 0<=frac<=1: raise ValueError("invalid descriptive WGBS coverage/fraction")
        for p in indexed_hits(promoter_index,promoter_max_length,r["chromosome"],s,e):
            if overlap(s,e,p["start"],p["end"]):
                sid=r.get("sample_id",NA); ov=overlap(s,e,p["start"],p["end"])
                cpg_links.append({"cpg_id":r["cpg_id"],"promoter_id":p["promoter_id"],"sample_id":sid,"overlap_bp":ov,"methylated_count":NA,"unmethylated_count":NA,"coverage":cov,"methylation_fraction":frac,"evidence_type":"DESCRIPTIVE_PERCENTAGE_WITH_COVERAGE"})
                z=meth_desc[(p["promoter_id"],sid)]; z[0]+=frac*cov; z[1]+=cov; z[2].add(r["cpg_id"])
    write_tsv(out/"relationships"/"cpg_promoter.tsv",cpg_links,["cpg_id","promoter_id","sample_id","overlap_bp","methylated_count","unmethylated_count","coverage","methylation_fraction","evidence_type"])
    meth_rows=[]
    for (pid,sid),(m,u,ids,npass) in sorted(meth.items()):
        cov=m+u; meth_rows.append({"promoter_id":pid,"sample_id":sid,"total_methylated_reads":m,"total_unmethylated_reads":u,"total_coverage":cov,"coverage_weighted_methylation_fraction":(m/cov if cov else NA),"number_of_observed_CpGs":len(ids),"number_of_CpGs_passing_analysis_coverage":npass,"evidence_type":"COUNT_PRESERVING"})
    for (pid,sid),(weighted,cov,ids) in sorted(meth_desc.items()):
        meth_rows.append({"promoter_id":pid,"sample_id":sid,"total_methylated_reads":NA,"total_unmethylated_reads":NA,"total_coverage":cov,"coverage_weighted_methylation_fraction":weighted/cov if cov else NA,"number_of_observed_CpGs":len(ids),"number_of_CpGs_passing_analysis_coverage":NA,"evidence_type":"DESCRIPTIVE_PERCENTAGE_WITH_COVERAGE"})
    write_tsv(out/"evidence"/"promoter_methylation_by_sample.tsv",meth_rows,["promoter_id","sample_id","total_methylated_reads","total_unmethylated_reads","total_coverage","coverage_weighted_methylation_fraction","number_of_observed_CpGs","number_of_CpGs_passing_analysis_coverage","evidence_type"])
    dmrs=loaded.get(("WGBS","wgbs_dmr"),[]); dmr_links=[]
    dmr_index,dmr_max_length=interval_index(dmrs) if dmrs else ({},0)
    for r in dmrs:
        s,e=int(r["start"]),int(r["end"])
        for p in indexed_hits(promoter_index,promoter_max_length,r["chromosome"],s,e):
            if overlap(s,e,p["start"],p["end"]):
                dmr_links.append({"dmr_id":r["dmr_id"],"promoter_id":p["promoter_id"],"gene_id_base":p["gene_id_base"],"overlap_bp":overlap(s,e,p["start"],p["end"]),"n_cpg":r.get("n_cpg",NA),"delta_methylation":r.get("delta_methylation",NA),"delta_methylation_pp":r.get("delta_methylation_pp",NA),"pvalue":r.get("pvalue",NA),"qvalue":r.get("qvalue",NA),"contrast":r.get("contrast",NA)})
    write_tsv(out/"relationships"/"dmr_promoter.tsv",dmr_links,["dmr_id","promoter_id","gene_id_base","overlap_bp","n_cpg","delta_methylation","delta_methylation_pp","pvalue","qvalue","contrast"])
    # ATAC/H3K27ac colocalization.
    coloc=[]
    adiff={r["region_id"]:r for r in loaded.get(("ATAC","atac_differential"),[])}; hdiff={r["region_id"]:r for r in loaded.get(("H3K27ac","h3k27ac_differential"),[])}
    chip_index,chip_max_length=interval_index(chip_regions) if chip_regions else ({},0)
    for ar in atac_regions:
        for hr in indexed_hits(chip_index,chip_max_length,ar["chromosome"],int(ar["start"]),int(ar["end"])):
            s=max(int(ar["start"]),int(hr["start"])); e=min(int(ar["end"]),int(hr["end"]))
            if e>s: coloc.append({"atac_region_id":ar["region_id"],"h3k27ac_region_id":hr["region_id"],"chromosome":ar["chromosome"],"overlap_start":s,"overlap_end":e,"overlap_bp":e-s,"atac_effect":adiff.get(ar["region_id"],{}).get("log2FoldChange",NA),"atac_status":state(adiff.get(ar["region_id"],{}).get("log2FoldChange"),adiff.get(ar["region_id"],{}).get("padj"),diff_complete.get("ATAC",False),a.fdr),"h3k27ac_effect":hdiff.get(hr["region_id"],{}).get("log2FoldChange",NA),"h3k27ac_status":state(hdiff.get(hr["region_id"],{}).get("log2FoldChange"),hdiff.get(hr["region_id"],{}).get("padj"),diff_complete.get("H3K27ac",False),a.fdr),"relationship":"ATAC_H3K27AC_COLOCATED"})
    write_tsv(out/"relationships"/"atac_h3k27ac_colocalization.tsv",coloc,["atac_region_id","h3k27ac_region_id","chromosome","overlap_start","overlap_end","overlap_bp","atac_effect","atac_status","h3k27ac_effect","h3k27ac_status","relationship"])
    # Summary table.
    by_atac=defaultdict(list); by_chip=defaultdict(list); by_dmr=defaultdict(list); by_meth=defaultdict(list)
    for x in atac_links: by_atac[x["promoter_id"]].append(x)
    for x in chip_links: by_chip[x["promoter_id"]].append(x)
    for x in dmr_links: by_dmr[x["promoter_id"]].append(x)
    for x in meth_rows: by_meth[x["promoter_id"]].append(x)
    summary=[]
    for p in promoters:
        gid=p["gene_id_base"]; rr=rna_diff.get(gid,{}); rb=rna_base.get(gid,{})
        rs=state(rr.get("log2FoldChange"),rr.get("padj"),diff_complete.get("RNA",False),a.fdr) if available["RNA"] else "NOT_AVAILABLE"
        al=by_atac[p["promoter_id"]]; hl=by_chip[p["promoter_id"]]; dl=by_dmr[p["promoter_id"]]
        astates=[state(adiff.get(x["source_region_id"],{}).get("log2FoldChange"),adiff.get(x["source_region_id"],{}).get("padj"),diff_complete.get("ATAC",False),a.fdr) for x in al]
        hstates=[state(hdiff.get(x["source_region_id"],{}).get("log2FoldChange"),hdiff.get(x["source_region_id"],{}).get("padj"),diff_complete.get("H3K27ac",False),a.fdr) for x in hl]
        dstates=[state(x.get("delta_methylation"),x.get("qvalue"),diff_complete.get("WGBS",False),a.fdr,True) for x in dl]
        def aggregate(states,available_flag,up="UP",down="DOWN"):
            if not available_flag:return "NOT_AVAILABLE"
            if not differential_allowed:return "NOT_TESTED"
            if up in states and down in states:return "MIXED_SIGNIFICANT"
            if up in states:return up
            if down in states:return down
            return "NO_SIGNIFICANT_CHANGE" if states else "NOT_TESTED"
        ats=aggregate(astates,available["ATAC"]); hts=aggregate(hstates,available["H3K27ac"]); mts=aggregate(dstates,available["WGBS"],"HYPERMETHYLATED","HYPOMETHYLATED")
        effects_a=[abs(fnum(adiff.get(x["source_region_id"],{}).get("log2FoldChange")) or 0) for x in al]
        effects_h=[abs(fnum(hdiff.get(x["source_region_id"],{}).get("log2FoldChange")) or 0) for x in hl]
        effects_m=[abs(fnum(x.get("delta_methylation")) or 0) for x in dl]
        methvals=by_meth[p["promoter_id"]]; count_meth=[x for x in methvals if x["evidence_type"]=="COUNT_PRESERVING"]
        mt=sum(int(x["total_methylated_reads"]) for x in count_meth); ut=sum(int(x["total_unmethylated_reads"]) for x in count_meth); cov=mt+ut
        if cov: promoter_fraction=mt/cov
        else:
            desc_cov=sum(int(x["total_coverage"]) for x in methvals)
            promoter_fraction=(sum(float(x["coverage_weighted_methylation_fraction"])*int(x["total_coverage"]) for x in methvals)/desc_cov) if desc_cov else NA
        pattern=f"RNA_{rs}__ATAC_{ats}__H3K27AC_{hts}__METHYLATION_{mts}"
        available_count=sum(available.values()); tested_count=sum(s not in {"NOT_AVAILABLE","NOT_TESTED"} for s in (rs,ats,hts,mts))
        contrast_text=",".join(next(iter(triples))) if len(triples)==1 else NA
        summary.append({**{k:p[k] for k in ("promoter_id","tss_id","gene_id_original","gene_id_base","gene_name","chromosome","promoter_start","promoter_end","strand")},"integration_match_type":a.match_type,"rna_available":"AVAILABLE" if available["RNA"] else "NOT_AVAILABLE","rna_effect":rr.get("log2FoldChange",NA),"rna_pvalue":rr.get("pvalue",NA),"rna_padj":rr.get("padj",NA),"rna_state":rs,"rna_descriptive_expression":rb.get("normalized_expression",rb.get("baseMean",NA)),"atac_available":"AVAILABLE" if available["ATAC"] else "NOT_AVAILABLE","n_atac_regions":len(al),"n_atac_up":astates.count("UP"),"n_atac_down":astates.count("DOWN"),"max_abs_atac_effect":max(effects_a) if effects_a else NA,"atac_state":ats,"wgbs_available":"AVAILABLE" if available["WGBS"] else "NOT_AVAILABLE","promoter_methylation_fraction":promoter_fraction,"n_observed_cpg":sum(int(x["number_of_observed_CpGs"]) for x in methvals),"n_dmr":len(dl),"n_hypermethylated_dmr":dstates.count("HYPERMETHYLATED"),"n_hypomethylated_dmr":dstates.count("HYPOMETHYLATED"),"max_abs_delta_methylation":max(effects_m) if effects_m else NA,"methylation_state":mts,"h3k27ac_available":"AVAILABLE" if available["H3K27ac"] else "NOT_AVAILABLE","n_h3k27ac_regions":len(hl),"n_h3k27ac_up":hstates.count("UP"),"n_h3k27ac_down":hstates.count("DOWN"),"max_abs_h3k27ac_effect":max(effects_h) if effects_h else NA,"h3k27ac_state":hts,"regulatory_pattern":pattern,"available_modality_count":available_count,"tested_modality_count":tested_count,"contrast":contrast_text,"integration_status":"COMPLETE","provenance_id":"INTEGRATION_PROVENANCE_001"})
    summary_fields=list(summary[0]) if summary else []
    write_tsv(out/"evidence"/"promoter_evidence.tsv",summary,summary_fields)
    # Gene summary preserves promoter heterogeneity.
    gs=[]
    for gid in sorted({r["gene_id_base"] for r in summary}):
        rr=[r for r in summary if r["gene_id_base"]==gid]; pats={r["regulatory_pattern"] for r in rr}
        gs.append({"gene_id_base":gid,"gene_id_original":safe_join(r["gene_id_original"] for r in rr),"gene_name":safe_join(r["gene_name"] for r in rr),"number_of_promoters":len(rr),"number_with_ATAC_evidence":sum(int(r["n_atac_regions"])>0 for r in rr),"number_with_H3K27ac_evidence":sum(int(r["n_h3k27ac_regions"])>0 for r in rr),"number_with_WGBS_evidence":sum(int(r["n_observed_cpg"])>0 for r in rr),"promoter_pattern_diversity":len(pats),"promoter_interpretation":"PROMOTER_HETEROGENEITY" if len(pats)>1 else "CONSISTENT_PROMOTER_PATTERN","gene_level_rna_state":rr[0]["rna_state"],"promoter_ids":safe_join(r["promoter_id"] for r in rr)})
    write_tsv(out/"evidence"/"gene_summary.tsv",gs,list(gs[0]) if gs else [])
    # Distal candidates, anchored on non-promoter ATAC or H3K27ac-only regions.
    promoter_region_ids={x["source_region_id"] for x in atac_links}; promoter_chip_ids={x["source_region_id"] for x in chip_links}
    anchors=[]
    for r in atac_regions:
        if r["region_id"] not in promoter_region_ids: anchors.append((r,"ATAC_ANCHORED"))
    if not atac_regions:
        for r in chip_regions:
            if r["region_id"] not in promoter_chip_ids: anchors.append((r,"H3K27AC_ONLY"))
    def distal_for(maxd):
        res=[]
        for r,anchor in anchors:
            s,e=int(r["start"]),int(r["end"]); candidates=[]
            for t in tss:
                if t["chromosome"]!=r["chromosome"]: continue
                tp=int(t["tss_start"]); d=0 if s<=tp<e else min(abs(tp-s),abs(tp-(e-1)))
                candidates.append((d,t))
            if not candidates: continue
            md=min(x[0] for x in candidates)
            if md>maxd: continue
            hids=[x["h3k27ac_region_id"] for x in coloc if x["atac_region_id"]==r["region_id"]]
            for d,t in candidates:
                if d==md: res.append({"region_id":r["region_id"],"candidate_gene_id":t["gene_id_base"],"candidate_tss_id":t["tss_id"],"distance_to_tss":d,"link_method":"NEAREST_TSS_CANDIDATE","max_distance":maxd,"source_modality":"ATAC" if anchor=="ATAC_ANCHORED" else "H3K27ac","anchor_type":anchor,"h3k27ac_region_ids":safe_join(hids)})
        return res
    distal=distal_for(a.distal_max)
    write_tsv(out/"relationships"/"distal_candidate_links.tsv",distal,["region_id","candidate_gene_id","candidate_tss_id","distance_to_tss","link_method","max_distance","source_modality","anchor_type","h3k27ac_region_ids"])
    # External maps remain separate.
    external=[]
    ext_path=a.external_map or next((r["path"] for r in manifest if r["data_type"]=="enhancer_gene_map"),None)
    if ext_path:
        for r in read_tsv(ext_path):
            if r.get("genome_build") not in {"GRCh38","GRCh38_primary","GRCh38.p14_primary"}: raise ValueError("external map genome build is incompatible")
            method="EXTERNAL_EXPERIMENTALLY_VALIDATED" if r.get("validation_type")=="EXPERIMENTALLY_VALIDATED" else "EXTERNAL_PREDICTED"
            external.append({**r,"link_method":method,"gene_id_base":base_gene(r.get("gene_identifier"))})
    ext_fields=["source","version","genome_build","chromosome","start","end","gene_identifier","gene_id_base","evidence_type","validation_type","link_method"]
    write_tsv(out/"relationships"/"external_enhancer_gene_links.tsv",external,ext_fields)
    # Sensitivity counts use the same strict overlap semantics; patterns are not optimized.
    sens=[]
    promoter_sets=[("-2000/+500",promoters)]+[(Path(x).stem.split("promoters.")[-1],read_tsv(x)) for x in a.sensitivity_promoters]
    primary_patterns={x["tss_id"]:x["regulatory_pattern"] for x in summary}
    def sensitivity_pattern(p):
        gid=p["gene_id_base"]; rr=rna_diff.get(gid,{})
        rs=state(rr.get("log2FoldChange"),rr.get("padj"),diff_complete.get("RNA",False),a.fdr) if available["RNA"] else "NOT_AVAILABLE"
        ars=indexed_hits(atac_index,atac_max_length,p["chromosome"],p["start"],p["end"])
        hrs=indexed_hits(chip_index,chip_max_length,p["chromosome"],p["start"],p["end"])
        drs=indexed_hits(dmr_index,dmr_max_length,p["chromosome"],p["start"],p["end"])
        av=[state(adiff.get(r["region_id"],{}).get("log2FoldChange"),adiff.get(r["region_id"],{}).get("padj"),diff_complete.get("ATAC",False),a.fdr) for r in ars]
        hv=[state(hdiff.get(r["region_id"],{}).get("log2FoldChange"),hdiff.get(r["region_id"],{}).get("padj"),diff_complete.get("H3K27ac",False),a.fdr) for r in hrs]
        mv=[state(r.get("delta_methylation"),r.get("qvalue"),diff_complete.get("WGBS",False),a.fdr,True) for r in drs]
        def agg(vals,present,positive,negative):
            if not present: return "NOT_AVAILABLE"
            if not differential_allowed: return "NOT_TESTED"
            if positive in vals and negative in vals: return "MIXED_SIGNIFICANT"
            if positive in vals: return positive
            if negative in vals: return negative
            return "NO_SIGNIFICANT_CHANGE" if vals else "NOT_TESTED"
        return f"RNA_{rs}__ATAC_{agg(av,available['ATAC'],'UP','DOWN')}__H3K27AC_{agg(hv,available['H3K27ac'],'UP','DOWN')}__METHYLATION_{agg(mv,available['WGBS'],'HYPERMETHYLATED','HYPOMETHYLATED')}"
    for label,ps in promoter_sets:
        for p in ps: p["start"]=int(p["promoter_start"]); p["end"]=int(p["promoter_end"])
        pidx,pmax=interval_index(ps)
        na=sum(len(indexed_hits(pidx,pmax,r["chromosome"],int(r["start"]),int(r["end"]))) for r in atac_regions)
        nh=sum(len(indexed_hits(pidx,pmax,r["chromosome"],int(r["start"]),int(r["end"]))) for r in chip_regions)
        alt={p["tss_id"]:sensitivity_pattern(p) for p in ps}
        changed=sum(alt.get(k)!=v for k,v in primary_patterns.items())
        sens.append({"analysis":"PROMOTER_WINDOW","setting":label,"atac_promoter_associations":na,"h3k27ac_promoter_associations":nh,"pattern_assignments_changed":changed,"distal_candidate_links":NA})
    for d in (50000,100000,250000): sens.append({"analysis":"DISTAL_MAX_DISTANCE","setting":str(d),"atac_promoter_associations":NA,"h3k27ac_promoter_associations":NA,"pattern_assignments_changed":NA,"distal_candidate_links":len(distal_for(d))})
    write_tsv(out/"sensitivity"/"integration.sensitivity_summary.tsv",sens,["analysis","setting","atac_promoter_associations","h3k27ac_promoter_associations","pattern_assignments_changed","distal_candidate_links"])
    # Machine-readable plot sources and deterministic SVG frequency chart.
    freq=defaultdict(int)
    for r in summary: freq[r["regulatory_pattern"]]+=1
    freqrows=[{"regulatory_pattern":k,"promoter_count":v} for k,v in sorted(freq.items())]
    write_tsv(out/"plots"/"regulatory_pattern_frequency.tsv",freqrows,["regulatory_pattern","promoter_count"])
    width=900; barh=22; height=50+barh*max(1,len(freqrows)); mx=max([r["promoter_count"] for r in freqrows] or [1])
    svg=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">','<rect width="100%" height="100%" fill="white"/>','<text x="10" y="22" font-family="sans-serif" font-size="16">Regulatory evidence-pattern frequency</text>']
    for i,r in enumerate(freqrows):
        y=40+i*barh; bw=int(300*r["promoter_count"]/mx); label=r["regulatory_pattern"]
        svg += [f'<rect x="580" y="{y-13}" width="{bw}" height="15" fill="#5B4BB7"/>',f'<text x="10" y="{y}" font-family="sans-serif" font-size="11">{label}</text>',f'<text x="{590+bw}" y="{y}" font-family="sans-serif" font-size="11">{r["promoter_count"]}</text>']
    svg.append('</svg>'); (out/"plots"/"regulatory_pattern_frequency.svg").write_text("\n".join(svg)+"\n",encoding="utf-8")
    if differential_allowed: differential_status="COMPLETE"
    elif a.mode=="baseline": differential_status="NOT_RUN_DISABLED"
    elif cstate=="CONTRAST_MISMATCH": differential_status="NOT_RUN_INVALID_DESIGN"
    else: differential_status="NOT_RUN_NO_CONTRAST"
    statuses=[{"component":"integration_reference","status":"COMPLETE"},{"component":"integration_promoter","status":"COMPLETE" if summary else "COMPLETE_EMPTY_RESULT"},{"component":"integration_differential","status":differential_status},{"component":"integration_distal","status":"COMPLETE" if distal else "COMPLETE_EMPTY_RESULT"}]
    write_tsv(out/"provenance"/"integration.analysis_status.tsv",statuses,["component","status"])
    provenance={"author":"Sarvesh Galgale","project":"regulatory-multiomics-pipeline","project_version":"0.1.0","reference_policy":{"species":"Homo sapiens","assembly":"GRCh38 primary assembly","annotation":"GENCODE 50","allowed_chromosomes":sorted(ALLOWED),"coordinate_system":"0-based-half-open"},"tss_derivation":"all distinct transcript TSS; exact duplicates collapsed only within gene and strand","promoter_primary":{"upstream":2000,"downstream":500},"promoter_sensitivity":[[-1000,500],[-2000,1000]],"distal_maximum_bp":a.distal_max,"distal_sensitivity_bp":[50000,100000,250000],"gene_id_policy":"terminal numeric Ensembl version removed for base join; original retained","overlap_semantics":"strict positive base-pair intersection; book-ended intervals do not overlap","significance_fdr":a.fdr,"integration_match_type":a.match_type,"contrast_compatibility":cstate,"no_cross_assay_pvalue_combination":True,"no_combined_regulatory_score":True,"inputs":checksums,"statuses":statuses}
    (out/"provenance"/"integration.provenance.json").write_text(json.dumps(provenance,sort_keys=True,indent=2)+"\n",encoding="utf-8")
if __name__=="__main__": main()
