#!/usr/bin/env python3
import argparse, gzip, hashlib, json, re
from pathlib import Path

PRIMARY=re.compile(r'^chr(?:[1-9]|1[0-9]|2[0-2]|X|Y|M)$')
def op(p): return gzip.open(p,'rt',encoding='utf-8') if str(p).endswith('.gz') else open(p,encoding='utf-8')
def digest(p,algo):
    h=hashlib.new(algo)
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()
def main():
    p=argparse.ArgumentParser(); p.add_argument('--fasta',required=True); p.add_argument('--gtf',required=True); p.add_argument('--scope',required=True); p.add_argument('--output',required=True); a=p.parse_args()
    f=Path(a.fasta); g=Path(a.gtf); contigs=[]; gtf=set(); comments=[]
    with op(f) as h:
        for line in h:
            if line.startswith('>'): contigs.append(line[1:].split()[0])
    with op(g) as h:
        for line in h:
            if line.startswith('#'):
                if len(comments)<30: comments.append(line.strip())
            elif line.strip(): gtf.add(line.split('\t',1)[0])
    errors=[]
    if not contigs: errors.append('FASTA contains no sequences')
    if not gtf: errors.append('GTF contains no features')
    if len(contigs)!=len(set(contigs)): errors.append('duplicate FASTA sequence names')
    if not gtf.issubset(set(contigs)): errors.append('GTF sequence names absent from FASTA: '+','.join(sorted(gtf-set(contigs))))
    nonprimary=sorted(x for x in set(contigs)|gtf if not PRIMARY.fullmatch(x))
    if nonprimary: errors.append('non-primary sequence names are not allowed: '+','.join(nonprimary[:10]))
    if a.scope=='production':
        expected={f'chr{i}' for i in range(1,23)}|{'chrX','chrY','chrM'}
        if set(contigs)!=expected: errors.append('production FASTA must contain exactly chr1-22, chrX, chrY, and chrM')
        if not any('version 50' in x.lower() or 'gencode v50' in x.lower() for x in comments): errors.append('production GTF does not declare GENCODE release 50')
    if a.scope=='compact_chr22_execution' and set(contigs)!={'chr22'}: errors.append('compact chr22 reference must contain only chr22')
    if errors: raise SystemExit('Reference validation failed:\n- '+'\n- '.join(errors))
    out={'status':'COMPLETE','reference_scope':a.scope,'coordinate_representation':'GRCh38 primary-assembly chromosome names','annotation_release':('GENCODE 50' if a.scope in {'production','compact_chr22_execution'} else 'synthetic'), 'fasta':{'name':f.name,'sha256':digest(f,'sha256'),'md5':digest(f,'md5'),'contigs':contigs},'gtf':{'name':g.name,'sha256':digest(g,'sha256'),'md5':digest(g,'md5'),'sequence_names':sorted(gtf)}}
    Path(a.output).write_text(json.dumps(out,sort_keys=True,indent=2)+'\n',encoding='utf-8')
if __name__=='__main__': main()
