#!/usr/bin/env python3
import argparse, csv, hashlib
from collections import defaultdict
from pathlib import Path
from wgbs_cpg_contract import fasta, cpg_start

NA={'','NA','N/A','.'}
def integer(value,name,row):
    try: result=int(value)
    except ValueError: raise SystemExit(f'row {row}: {name} must be an integer')
    if result<0: raise SystemExit(f'row {row}: {name} must be nonnegative')
    return result

def main():
    p=argparse.ArgumentParser(); p.add_argument('--input',required=True); p.add_argument('--fasta',required=True); p.add_argument('--sample-id',required=True); p.add_argument('--library-id',required=True); p.add_argument('--type',required=True,choices=['cpg_count_tsv','cpg_percentage_tsv']); p.add_argument('--output',required=True); a=p.parse_args()
    seq=fasta(a.fasta); merged=defaultdict(lambda:[0,0,[],set()]); source_coordinates=set()
    with open(a.input,encoding='utf-8-sig',newline='') as handle:
        reader=csv.DictReader(handle,delimiter='\t')
        required={'chromosome','start','end','coordinate_system','methylation_fraction'}
        if a.type=='cpg_count_tsv': required|={'methylated_count'}
        if not reader.fieldnames or required-set(reader.fieldnames): raise SystemExit('processed WGBS file is missing required columns: '+', '.join(sorted(required-set(reader.fieldnames or []))))
        for row_no,row in enumerate(reader,2):
            chrom=row['chromosome']; coord=row['coordinate_system']; start=integer(row['start'],'start',row_no); end=integer(row['end'],'end',row_no)
            source_coordinates.add(coord)
            if coord=='1-based_closed_cytosine':
                if end != start: raise SystemExit(f'row {row_no}: 1-based cytosine interval must have identical start and end')
                source=start-1
            elif coord=='0-based-half-open':
                if end <= start: raise SystemExit(f'row {row_no}: end must be greater than start')
                if end-start not in {1,2}: raise SystemExit(f'row {row_no}: 0-based CpG/cytosine interval must span one or two bases')
                source=start
            else: raise SystemExit(f'row {row_no}: unsupported coordinate_system {coord}')
            canonical,strand=cpg_start(seq,chrom,source)
            frac=row.get('methylation_fraction','')
            if frac not in NA:
                try: frac=float(frac)
                except ValueError: raise SystemExit(f'row {row_no}: methylation_fraction must be numeric or NA')
                if not 0<=frac<=1: raise SystemExit(f'row {row_no}: methylation_fraction outside [0,1]')
            if a.type=='cpg_count_tsv':
                meth=integer(row['methylated_count'],'methylated_count',row_no)
                if 'unmethylated_count' in row and row['unmethylated_count'] not in NA:
                    unmeth=integer(row['unmethylated_count'],'unmethylated_count',row_no)
                    if 'coverage' in row and row['coverage'] not in NA and integer(row['coverage'],'coverage',row_no)!=meth+unmeth: raise SystemExit(f'row {row_no}: coverage conflicts with methylated and unmethylated counts')
                elif 'coverage' in row and row['coverage'] not in NA:
                    coverage=integer(row['coverage'],'coverage',row_no)
                    if meth>coverage: raise SystemExit(f'row {row_no}: methylated_count exceeds coverage')
                    unmeth=coverage-meth
                else: raise SystemExit(f'row {row_no}: count input needs unmethylated_count or coverage')
                coverage=meth+unmeth
                if coverage==0 and frac not in NA: raise SystemExit(f'row {row_no}: zero coverage must have NA methylation_fraction')
                if coverage and frac not in NA and abs(frac-meth/coverage)>1e-6: raise SystemExit(f'row {row_no}: methylation_fraction conflicts with counts')
                merged[(chrom,canonical,canonical+2)][0]+=meth; merged[(chrom,canonical,canonical+2)][1]+=unmeth
            else:
                if frac in NA: merged[(chrom,canonical,canonical+2)][2].append(None)
                else: merged[(chrom,canonical,canonical+2)][2].append(frac)
            merged[(chrom,canonical,canonical+2)][3].add(strand)
    source_coordinate_system=','.join(sorted(source_coordinates))
    fields=['cpg_id','chromosome','start','end','coordinate_system','sample_id','library_id','methylated_count','unmethylated_count','coverage','methylation_fraction','source_strands','measurement_state','source_coordinate_system']
    with open(a.output,'w',newline='',encoding='utf-8') as out:
        writer=csv.DictWriter(out,fieldnames=fields,delimiter='\t',lineterminator='\n'); writer.writeheader()
        for (chrom,start,end),(meth,unmeth,fractions,strands) in sorted(merged.items()):
            cid='CpG_'+hashlib.sha256(f'{chrom}:{start}-{end}'.encode()).hexdigest()[:16]
            if a.type=='cpg_count_tsv': cov=meth+unmeth; fraction=f'{meth/cov:.8f}' if cov else 'NA'; state='OBSERVED' if cov else 'ZERO_COVERAGE'
            else: cov='NA'; meth='NA'; unmeth='NA'; valid=[x for x in fractions if x is not None]; fraction=f'{sum(valid)/len(valid):.8f}' if valid else 'NA'; state='DESCRIPTIVE_PERCENTAGE_ONLY' if valid else 'MISSING'
            writer.writerow({'cpg_id':cid,'chromosome':chrom,'start':start,'end':end,'coordinate_system':'0-based-half-open','sample_id':a.sample_id,'library_id':a.library_id,'methylated_count':meth,'unmethylated_count':unmeth,'coverage':cov,'methylation_fraction':fraction,'source_strands':','.join(sorted(strands)),'measurement_state':state,'source_coordinate_system':'processed_'+source_coordinate_system})

if __name__=='__main__': main()
