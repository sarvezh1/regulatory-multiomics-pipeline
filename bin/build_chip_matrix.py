#!/usr/bin/env python3
import argparse,csv
p=argparse.ArgumentParser(); p.add_argument('--regions',required=True); p.add_argument('--evidence',nargs='+',required=True); p.add_argument('--output',required=True); a=p.parse_args()
regions=list(csv.DictReader(open(a.regions),delimiter='\t')); samples=[]; data=[]
for path in a.evidence:
    sid=path.rsplit('/',1)[-1].split('.evidence.bed')[0]; samples.append(sid); d={}
    for line in open(path):
        x=line.rstrip().split('\t'); d.setdefault(x[0],[]).append((int(x[1]),int(x[2])))
    data.append(d)
with open(a.output,'w',newline='') as out:
    w=csv.writer(out,delimiter='\t',lineterminator='\n'); w.writerow(['region_id','chromosome','start','end']+samples)
    for r in regions:
        s,e=int(r['start']),int(r['end']); w.writerow([r['region_id'],r['chromosome'],s,e]+[sum(x<e and y>s for x,y in d.get(r['chromosome'],[])) for d in data])
