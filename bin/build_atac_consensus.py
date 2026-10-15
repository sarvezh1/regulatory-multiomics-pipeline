#!/usr/bin/env python3
import argparse, hashlib
def main():
    p=argparse.ArgumentParser(); p.add_argument('--peaks',nargs='+',required=True); p.add_argument('--regions',required=True); p.add_argument('--provenance',required=True); a=p.parse_args()
    peaks=[]
    for path in a.peaks:
        sample=path.split('/')[-1].split('.')[0]
        with open(path) as h:
            for n,line in enumerate(h,1):
                if not line.strip() or line.startswith('#'): continue
                x=line.rstrip().split('\t'); s,e=int(x[1]),int(x[2])
                if s<0 or e<=s: raise SystemExit(f'invalid peak coordinates in {path}:{n}')
                peaks.append((x[0],s,e,sample,x[3] if len(x)>3 else f'{sample}_peak_{n}'))
    peaks.sort(key=lambda z:(z[0],z[1],z[2],z[3],z[4])); merged=[]
    for peak in peaks:
        if not merged or peak[0]!=merged[-1][0] or peak[1]>=merged[-1][2]: merged.append([peak[0],peak[1],peak[2],[peak]])
        else: merged[-1][2]=max(merged[-1][2],peak[2]); merged[-1][3].append(peak)
    with open(a.regions,'w') as out, open(a.provenance,'w') as prov:
        out.write('region_id\tchromosome\tstart\tend\tcoordinate_system\n'); prov.write('region_id\tsample_id\tsource_peak_id\tsource_chromosome\tsource_start\tsource_end\n')
        for chrom,s,e,contributors in merged:
            rid='ATAC_'+hashlib.sha256(f'{chrom}:{s}-{e}'.encode()).hexdigest()[:16]
            out.write(f'{rid}\t{chrom}\t{s}\t{e}\t0-based-half-open\n')
            for p in contributors: prov.write(f'{rid}\t{p[3]}\t{p[4]}\t{p[0]}\t{p[1]}\t{p[2]}\n')
if __name__=='__main__': main()
