#!/usr/bin/env python3
import argparse,csv,json
p=argparse.ArgumentParser(); p.add_argument('--sample',required=True); p.add_argument('--metrics',nargs='+',required=True); p.add_argument('--peaks',required=True); p.add_argument('--evidence',required=True); p.add_argument('--output',required=True); a=p.parse_args()
vals={}
for path in a.metrics:
    for r in csv.DictReader(open(path),delimiter='\t'):
        try: vals[r['metric']]=vals.get(r['metric'],0)+int(r['value'])
        except ValueError: vals[r['metric']]=r['value']
peaks=sum(1 for x in open(a.peaks) if x.strip() and not x.startswith('#'))
ev=[]
for line in open(a.evidence):
    x=line.rstrip().split('\t'); ev.append((x[0],int(x[1]),int(x[2])))
pk=[]
for line in open(a.peaks):
    if line.strip() and not line.startswith('#'):
        x=line.rstrip().split('\t'); pk.append((x[0],int(x[1]),int(x[2])))
in_peak=sum(any(c==pc and s<pe and e>ps for pc,ps,pe in pk) for c,s,e in ev)
usable=len(ev); vals.update({'sample_id':a.sample,'peak_count':peaks,'usable_evidence':usable,'reads_or_fragments_in_peaks':in_peak,'FRiP':in_peak/usable if usable else None,'FRiP_denominator':'usable nonduplicate nuclear reads/fragments after MAPQ and blacklist filtering','NSC':'NOT_IMPLEMENTED','RSC':'NOT_IMPLEMENTED'})
open(a.output,'w').write(json.dumps(vals,sort_keys=True,indent=2)+'\n')
