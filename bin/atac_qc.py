#!/usr/bin/env python3
import argparse, csv, json, math
from collections import Counter

def main():
    p=argparse.ArgumentParser(); p.add_argument('--fragments',required=True); p.add_argument('--pre-blacklist-metrics',required=True); p.add_argument('--filter-metrics',required=True); p.add_argument('--peaks',required=True); p.add_argument('--tss',required=True); p.add_argument('--output',required=True); p.add_argument('--lengths',required=True); a=p.parse_args()
    fragments=[]
    with open(a.fragments) as h:
        for line in h:
            x=line.rstrip().split('\t'); fragments.append((x[0],int(x[1]),int(x[2])))
    peaks={}
    with open(a.peaks) as h:
        for line in h:
            if not line.strip(): continue
            x=line.rstrip().split('\t'); peaks.setdefault(x[0],[]).append((int(x[1]),int(x[2])))
    lengths=Counter(e-s for _,s,e in fragments)
    with open(a.lengths,'w') as out:
        out.write('fragment_length\tcount\n')
        for length,count in sorted(lengths.items()): out.write(f'{length}\t{count}\n')
    overlap=sum(any(s<pe and e>ps for ps,pe in peaks.get(c,[])) for c,s,e in fragments)
    insertions={}
    for c,s,e in fragments: insertions.setdefault(c,[]).extend((s+4,e-5))
    centers=[]; flanks=[]
    with open(a.tss) as h:
        for line in h:
            x=line.rstrip().split('\t'); c=x[0]; t=int(x[1]); vals=insertions.get(c,[])
            centers.append(sum(t-50 <= z < t+50 for z in vals)/100)
            flanks.append((sum(t-1000 <= z < t-900 for z in vals)+sum(t+900 <= z < t+1000 for z in vals))/200)
    center=sum(centers)/len(centers) if centers else 0; flank=sum(flanks)/len(flanks) if flanks else 0
    tss_enrichment=(center/flank) if flank else None
    pre=json.load(open(a.pre_blacklist_metrics))
    with open(a.filter_metrics) as h:
        filt={r['metric']:r['value'] for r in csv.DictReader(h,delimiter='\t')}
    def fraction(n,d): return (float(n)/float(d)) if float(d) else None
    total=len(fragments); nfr=sum((e-s)<150 for _,s,e in fragments); mono=sum(150<=(e-s)<300 for _,s,e in fragments); multi=sum((e-s)>=300 for _,s,e in fragments)
    mt_available=filt['mitochondrial_reference_available']=='1'
    qc={'total_fragments':int(filt['all_first_mates']),'mapped_fraction':fraction(filt['mapped_first_mates'],filt['all_first_mates']),'properly_paired_fraction':fraction(filt['proper_primary'],filt['all_first_mates']),'mitochondrial_reference_available':mt_available,'mitochondrial_fraction':fraction(filt['mitochondrial'],filt['mapped_first_mates']) if mt_available else None,'duplicate_fraction':fraction(filt['duplicate'],filt['mapped_first_mates']),'usable_fragments':total,'blacklist_fragments':pre['blacklist_fragments'],'frip':overlap/total if total else None,'fragments_in_peaks':overlap,'peak_count':sum(len(v) for v in peaks.values()),'nucleosome_free_fragments':nfr,'mononucleosomal_fragments':mono,'multinucleosomal_fragments':multi,'tss_enrichment':tss_enrichment,'tss_positions':len(centers),'NRF':None if filt['NRF']=='NA' else float(filt['NRF']),'PBC1':None if filt['PBC1']=='NA' else float(filt['PBC1']),'PBC2':None if filt['PBC2']=='NA' else float(filt['PBC2']),'coordinate_policy':'canonical BAM unchanged; insertion coordinates derived as +4/-5 for TSS QC only'}
    with open(a.output,'w') as out: json.dump(qc,out,sort_keys=True,indent=2); out.write('\n')
if __name__=='__main__': main()
