#!/usr/bin/env python3
import argparse, hashlib
p=argparse.ArgumentParser(); p.add_argument('--peaks',nargs='+',required=True); p.add_argument('--regions',required=True); p.add_argument('--provenance',required=True); a=p.parse_args()
peaks=[]
for path in a.peaks:
    sample=path.rsplit('/',1)[-1].split('_peaks.')[0]
    for n,line in enumerate(open(path),1):
        if not line.strip() or line.startswith('#'): continue
        x=line.rstrip().split('\t'); s,e=int(x[1]),int(x[2])
        if s<0 or e<=s: raise SystemExit(f'invalid peak coordinates in {path}:{n}')
        peaks.append((x[0],s,e,sample,x[3] if len(x)>3 else f'{sample}_peak_{n}'))
peaks.sort(); merged=[]
for x in peaks:
    if not merged or x[0]!=merged[-1][0] or x[1]>=merged[-1][2]: merged.append([x[0],x[1],x[2],[x]])
    else: merged[-1][2]=max(merged[-1][2],x[2]); merged[-1][3].append(x)
with open(a.regions,'w') as out, open(a.provenance,'w') as prov:
    out.write('region_id\tchromosome\tstart\tend\tcoordinate_system\n'); prov.write('region_id\tsample_id\tsource_peak_id\tsource_chromosome\tsource_start\tsource_end\n')
    for chrom,s,e,src in merged:
        rid='H3K27AC_'+hashlib.sha256(f'{chrom}:{s}-{e}'.encode()).hexdigest()[:16]
        out.write(f'{rid}\t{chrom}\t{s}\t{e}\t0-based-half-open\n')
        for x in src: prov.write(f'{rid}\t{x[3]}\t{x[4]}\t{x[0]}\t{x[1]}\t{x[2]}\n')
