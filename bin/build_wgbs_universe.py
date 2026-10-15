#!/usr/bin/env python3
import argparse, csv

def main():
    p=argparse.ArgumentParser(); p.add_argument('--inputs',nargs='+',required=True); p.add_argument('--metadata',required=True); p.add_argument('--min-coverage',type=int,default=5); p.add_argument('--output',required=True); p.add_argument('--summary',required=True); a=p.parse_args()
    with open(a.metadata,encoding='utf-8') as handle: samples=sorted({row['sample_id'] for row in csv.DictReader(handle,delimiter='\t')})
    data={}; loci={}
    for path in a.inputs:
        with open(path,encoding='utf-8') as handle:
            for row in csv.DictReader(handle,delimiter='\t'):
                sid=row['sample_id']; key=row['cpg_id']; loci[key]=(row['chromosome'],row['start'],row['end'],row['coordinate_system']); data[(key,sid)]=row
    fields=['cpg_id','chromosome','start','end','coordinate_system','sample_id','methylated_count','unmethylated_count','coverage','methylation_fraction','measurement_state','passes_min_coverage']
    before=len(loci); passing=set(); missing={s:0 for s in samples}
    with open(a.output,'w',newline='',encoding='utf-8') as out:
        writer=csv.DictWriter(out,fieldnames=fields,delimiter='\t',lineterminator='\n'); writer.writeheader()
        for cid,(chrom,start,end,coord) in sorted(loci.items(),key=lambda x:(x[1][0],int(x[1][1]))):
            for sid in samples:
                row=data.get((cid,sid))
                if row is None:
                    values={'methylated_count':'NA','unmethylated_count':'NA','coverage':'NA','methylation_fraction':'NA','measurement_state':'MISSING'}; missing[sid]+=1; passes='false'
                else:
                    values={k:row[k] for k in ('methylated_count','unmethylated_count','coverage','methylation_fraction','measurement_state')}
                    passes=str(row['coverage']!='NA' and int(row['coverage'])>=a.min_coverage).lower()
                    if passes=='true': passing.add(cid)
                writer.writerow({'cpg_id':cid,'chromosome':chrom,'start':start,'end':end,'coordinate_system':coord,'sample_id':sid,**values,'passes_min_coverage':passes})
    with open(a.summary,'w',newline='',encoding='utf-8') as out:
        writer=csv.writer(out,delimiter='\t',lineterminator='\n'); writer.writerow(['metric','value']); writer.writerow(['minimum_coverage',a.min_coverage]); writer.writerow(['cpg_loci_before_filtering',before]); writer.writerow(['cpg_loci_passing_in_at_least_one_sample',len(passing)])
        for sid in samples: writer.writerow([f'missing_loci_{sid}',missing[sid]])

if __name__=='__main__': main()
