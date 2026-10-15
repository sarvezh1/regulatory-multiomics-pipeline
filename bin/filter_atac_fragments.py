#!/usr/bin/env python3
import argparse, json
from pathlib import Path
def load_bed(path):
    d={}
    with open(path) as h:
        for n,line in enumerate(h,1):
            if not line.strip() or line.startswith('#'): continue
            x=line.rstrip().split('\t');
            if len(x)<3: raise SystemExit(f'malformed BED line {n}')
            s,e=int(x[1]),int(x[2])
            if s<0 or e<=s: raise SystemExit(f'invalid BED coordinates at line {n}')
            d.setdefault(x[0],[]).append((s,e))
    return d
def main():
    p=argparse.ArgumentParser(); p.add_argument('--fragments',required=True); p.add_argument('--blacklist',required=True); p.add_argument('--output',required=True); p.add_argument('--names',required=True); p.add_argument('--metrics',required=True); a=p.parse_args()
    bl=load_bed(a.blacklist); total=removed=0
    with open(a.fragments) as src, open(a.output,'w') as out, open(a.names,'w') as names:
        for n,line in enumerate(src,1):
            x=line.rstrip().split('\t')
            if len(x)<4: raise SystemExit(f'malformed fragment line {n}')
            s,e=int(x[1]),int(x[2]); total+=1
            if any(s < be and e > bs for bs,be in bl.get(x[0],[])): removed+=1; continue
            out.write(line); names.write(x[3]+'\n')
    Path(a.metrics).write_text(json.dumps({'pre_blacklist_fragments':total,'blacklist_fragments':removed,'usable_fragments':total-removed},sort_keys=True,indent=2)+'\n')
if __name__=='__main__': main()
