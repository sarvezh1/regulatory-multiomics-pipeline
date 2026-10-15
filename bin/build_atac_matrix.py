#!/usr/bin/env python3
import argparse, csv
def main():
    p=argparse.ArgumentParser(); p.add_argument('--regions',required=True); p.add_argument('--fragments',nargs='+',required=True); p.add_argument('--output',required=True); a=p.parse_args()
    with open(a.regions) as h: regions=list(csv.DictReader(h,delimiter='\t'))
    samples=[]; per=[]
    for path in a.fragments:
        sample=path.split('/')[-1].split('.')[0]; samples.append(sample); d={}
        with open(path) as h:
            for line in h:
                x=line.rstrip().split('\t'); d.setdefault(x[0],[]).append((int(x[1]),int(x[2])))
        per.append(d)
    with open(a.output,'w',newline='') as out:
        w=csv.writer(out,delimiter='\t',lineterminator='\n'); w.writerow(['region_id','chromosome','start','end']+samples)
        for r in regions:
            s,e=int(r['start']),int(r['end']); counts=[sum(fs<e and fe>s for fs,fe in d.get(r['chromosome'],[])) for d in per]
            w.writerow([r['region_id'],r['chromosome'],s,e]+counts)
if __name__=='__main__': main()
