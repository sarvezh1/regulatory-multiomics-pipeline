#!/usr/bin/env python3
import argparse, hashlib, re
from pathlib import Path

def attrs(text):
    return dict(re.findall(r'(\S+) "([^"]+)"',text))
def main():
    p=argparse.ArgumentParser(); p.add_argument('--gtf',required=True); p.add_argument('--bed',required=True); p.add_argument('--provenance',required=True); a=p.parse_args()
    positions={}
    with open(a.gtf,encoding='utf-8') as h:
        for line in h:
            if line.startswith('#'): continue
            x=line.rstrip('\n').split('\t')
            if len(x)!=9 or x[2]!='transcript': continue
            at=attrs(x[8]); pos=int(x[3])-1 if x[6]=='+' else int(x[4])-1
            key=(x[0],pos,x[6]); positions.setdefault(key,[]).append((at.get('gene_id','UNKNOWN'),at.get('transcript_id','UNKNOWN')))
    with open(a.bed,'w',encoding='utf-8') as out, open(a.provenance,'w',encoding='utf-8') as prov:
        prov.write('tss_id\tchromosome\tposition_0based\tstrand\tgene_ids\ttranscript_ids\n')
        for chrom,pos,strand in sorted(positions,key=lambda z:(z[0],z[1],z[2])):
            records=positions[(chrom,pos,strand)]; token=hashlib.sha1(f'{chrom}:{pos}:{strand}'.encode()).hexdigest()[:12]; tid=f'TSS_{token}'
            out.write(f'{chrom}\t{pos}\t{pos+1}\t{tid}\t0\t{strand}\n')
            prov.write(f'{tid}\t{chrom}\t{pos}\t{strand}\t{";".join(sorted({r[0] for r in records}))}\t{";".join(sorted({r[1] for r in records}))}\n')
if __name__=='__main__': main()
