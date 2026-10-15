#!/usr/bin/env python3
import argparse, csv, hashlib
import pyBigWig

p=argparse.ArgumentParser(); p.add_argument("--url",required=True); p.add_argument("--output",required=True); p.add_argument("--sample-id",default="ENCSR765JPC_rep1"); a=p.parse_args()
bw=pyBigWig.open(a.url)
if not bw or not bw.isBigBed(): raise SystemExit("remote WGBS resource is not a readable bigBed")
length=bw.chroms().get("chr22")
if not length: raise SystemExit("chr22 is absent from WGBS bigBed")
dyads={}
for left in range(0,length,1000000):
    entries=bw.entries("chr22",left,min(length,left+1000000)) or []
    for start,end,rest in entries:
        fields=rest.split("\t")
        strand=fields[2]; coverage=int(fields[6]); percent=float(fields[7]); context=fields[8]
        if context!="CG" or coverage<0 or not 0<=percent<=100: continue
        canonical=start if strand=="+" else start-1
        if canonical<0: continue
        z=dyads.setdefault(canonical,[0.0,0,set()]); z[0]+=coverage*(percent/100.0); z[1]+=coverage; z[2].add(strand)
bw.close()
fields=["cpg_id","chromosome","start","end","coordinate_system","sample_id","coverage","methylation_fraction","measurement_state","source_strands","exact_counts_available"]
with open(a.output,"w",encoding="utf-8",newline="") as h:
    w=csv.DictWriter(h,fieldnames=fields,delimiter="\t",lineterminator="\n"); w.writeheader()
    for pos,(weighted,cov,strands) in sorted(dyads.items()):
        token=hashlib.sha256(f"chr22:{pos}-{pos+2}".encode()).hexdigest()[:16]
        w.writerow({"cpg_id":"CpG_"+token,"chromosome":"chr22","start":pos,"end":pos+2,"coordinate_system":"0-based-half-open","sample_id":a.sample_id,"coverage":cov,"methylation_fraction":format(weighted/cov,".8f") if cov else "NA","measurement_state":"DESCRIPTIVE_PERCENTAGE_WITH_COVERAGE","source_strands":";".join(sorted(strands)),"exact_counts_available":"false"})
print(f"canonical_chr22_cpgs={len(dyads)}")
