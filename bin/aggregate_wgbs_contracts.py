#!/usr/bin/env python3
import argparse, csv
from collections import defaultdict

def main():
    p=argparse.ArgumentParser(); p.add_argument('--sample-id',required=True); p.add_argument('--inputs',nargs='+',required=True); p.add_argument('--output',required=True); a=p.parse_args()
    merged={}; libraries=set(); percentage_only=False
    for path in a.inputs:
        with open(path,encoding='utf-8') as handle:
            for row in csv.DictReader(handle,delimiter='\t'):
                libraries.add(row['library_id']); key=(row['cpg_id'],row['chromosome'],int(row['start']),int(row['end']))
                item=merged.setdefault(key,{'meth':0,'unmeth':0,'fractions':[],'strands':set(),'states':set()})
                item['states'].add(row['measurement_state']); item['strands'].update(filter(None,row['source_strands'].split(',')))
                if row['methylated_count']=='NA':
                    percentage_only=True
                    if row['methylation_fraction']!='NA': item['fractions'].append(float(row['methylation_fraction']))
                else: item['meth']+=int(row['methylated_count']); item['unmeth']+=int(row['unmethylated_count'])
    fields=['cpg_id','chromosome','start','end','coordinate_system','sample_id','contributing_libraries','methylated_count','unmethylated_count','coverage','methylation_fraction','source_strands','measurement_state']
    with open(a.output,'w',newline='',encoding='utf-8') as out:
        writer=csv.DictWriter(out,fieldnames=fields,delimiter='\t',lineterminator='\n'); writer.writeheader()
        for (cid,chrom,start,end),item in sorted(merged.items(),key=lambda x:(x[0][1],x[0][2],x[0][3])):
            if percentage_only:
                meth=unmeth=cov='NA'; fraction=f'{sum(item["fractions"])/len(item["fractions"]):.8f}' if item['fractions'] else 'NA'; state='DESCRIPTIVE_PERCENTAGE_ONLY' if item['fractions'] else 'MISSING'
            else:
                meth=item['meth']; unmeth=item['unmeth']; cov=meth+unmeth; fraction=f'{meth/cov:.8f}' if cov else 'NA'; state='OBSERVED' if cov else 'ZERO_COVERAGE'
            writer.writerow({'cpg_id':cid,'chromosome':chrom,'start':start,'end':end,'coordinate_system':'0-based-half-open','sample_id':a.sample_id,'contributing_libraries':','.join(sorted(libraries)),'methylated_count':meth,'unmethylated_count':unmeth,'coverage':cov,'methylation_fraction':fraction,'source_strands':','.join(sorted(item['strands'])),'measurement_state':state})

if __name__=='__main__': main()
