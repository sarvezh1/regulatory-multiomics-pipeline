#!/usr/bin/env python3
import argparse, csv, gzip, hashlib
from collections import defaultdict
from pathlib import Path

def open_text(path):
    return gzip.open(path,'rt') if str(path).endswith('.gz') else open(path,encoding='utf-8')

def fasta(path):
    sequences={}; name=None; parts=[]
    with open(path,encoding='utf-8') as handle:
        for line in handle:
            if line.startswith('>'):
                if name is not None: sequences[name]=''.join(parts).upper()
                name=line[1:].split()[0]; parts=[]
            else: parts.append(line.strip())
    if name is not None: sequences[name]=''.join(parts).upper()
    return sequences

def cpg_start(sequences,chrom,start):
    if chrom not in sequences or start < 0 or start >= len(sequences[chrom]): raise ValueError(f'coordinate outside reference: {chrom}:{start}')
    seq=sequences[chrom]
    if seq[start:start+2]=='CG': return start,'top'
    if start>0 and seq[start-1:start+1]=='CG': return start-1,'bottom'
    raise ValueError(f'non-CpG source coordinate: {chrom}:{start}')

def main():
    p=argparse.ArgumentParser(); p.add_argument('--coverage',required=True); p.add_argument('--fasta',required=True); p.add_argument('--sample-id',required=True); p.add_argument('--library-id',required=True); p.add_argument('--output',required=True); p.add_argument('--rejected-output',required=True); a=p.parse_args()
    seq=fasta(a.fasta); merged=defaultdict(lambda:[0,0,set()]); rejected=[]
    with open_text(a.coverage) as handle:
        for line_no,line in enumerate(handle,1):
            if not line.strip() or line.startswith('#'): continue
            fields=line.rstrip('\n').split('\t')
            if len(fields)<6: raise SystemExit(f'invalid Bismark coverage row {line_no}')
            chrom,start_s,end_s,percent,meth_s,unmeth_s=fields[:6]
            source_position=int(start_s); source_end=int(end_s); meth=int(meth_s); unmeth=int(unmeth_s)
            if source_position<1 or source_end!=source_position or meth<0 or unmeth<0: raise SystemExit(f'invalid 1-based Bismark coverage values at row {line_no}')
            try: canonical,strand=cpg_start(seq,chrom,source_position-1)
            except ValueError as error:
                rejected.append({'chromosome':chrom,'source_position_1based':source_position,'methylated_count':meth,'unmethylated_count':unmeth,'reason':str(error)})
                continue
            key=(chrom,canonical,canonical+2); merged[key][0]+=meth; merged[key][1]+=unmeth; merged[key][2].add(strand)
    header=['cpg_id','chromosome','start','end','coordinate_system','sample_id','library_id','methylated_count','unmethylated_count','coverage','methylation_fraction','source_strands','measurement_state','source_coordinate_system']
    with open(a.output,'w',newline='',encoding='utf-8') as out:
        writer=csv.DictWriter(out,fieldnames=header,delimiter='\t',lineterminator='\n'); writer.writeheader()
        for (chrom,start,end),(meth,unmeth,strands) in sorted(merged.items()):
            cov=meth+unmeth; key=f'{chrom}:{start}-{end}'; cid='CpG_'+hashlib.sha256(key.encode()).hexdigest()[:16]
            writer.writerow({'cpg_id':cid,'chromosome':chrom,'start':start,'end':end,'coordinate_system':'0-based-half-open','sample_id':a.sample_id,'library_id':a.library_id,'methylated_count':meth,'unmethylated_count':unmeth,'coverage':cov,'methylation_fraction':f'{meth/cov:.8f}' if cov else 'NA','source_strands':','.join(sorted(strands)),'measurement_state':'OBSERVED' if cov else 'ZERO_COVERAGE','source_coordinate_system':'Bismark_1_based_closed_coverage'})
    with open(a.rejected_output,'w',newline='',encoding='utf-8') as out:
        fields=['chromosome','source_position_1based','methylated_count','unmethylated_count','reason']; writer=csv.DictWriter(out,fieldnames=fields,delimiter='\t',lineterminator='\n'); writer.writeheader(); writer.writerows(rejected)

if __name__=='__main__': main()
