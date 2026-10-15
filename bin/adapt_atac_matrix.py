#!/usr/bin/env python3
import argparse, csv, hashlib, json
def main():
    p=argparse.ArgumentParser(); p.add_argument('--matrix',required=True); p.add_argument('--metadata',required=True); p.add_argument('--output',required=True); p.add_argument('--regions',required=True); p.add_argument('--provenance',required=True); a=p.parse_args()
    with open(a.metadata,encoding='utf-8') as h: meta=list(csv.DictReader(h,delimiter='\t'))
    sample_ids=[r['sample_id'] for r in meta]
    if len(sample_ids)!=len(set(sample_ids)): raise SystemExit('processed route requires one row per biological sample')
    with open(a.matrix,encoding='utf-8') as h: reader=csv.DictReader(h,delimiter='\t'); fields=reader.fieldnames or []; rows=list(reader)
    need=['region_id','chromosome','start','end']
    if any(x not in fields for x in need): raise SystemExit('region matrix requires region_id, chromosome, start, end')
    if any(x not in fields for x in sample_ids): raise SystemExit('region matrix lacks one or more samplesheet sample columns')
    seen=set()
    with open(a.output,'w',newline='') as out, open(a.regions,'w',newline='') as reg:
        wo=csv.writer(out,delimiter='\t',lineterminator='\n'); wr=csv.writer(reg,delimiter='\t',lineterminator='\n')
        wo.writerow(need+sample_ids); wr.writerow(need+['coordinate_system'])
        for n,r in enumerate(rows,2):
            try: s,e=int(r['start']),int(r['end'])
            except: raise SystemExit(f'non-integer coordinates at line {n}')
            if s<0 or e<=s: raise SystemExit(f'invalid coordinates at line {n}')
            if r['region_id'] in seen: raise SystemExit(f'duplicate region_id {r["region_id"]}')
            seen.add(r['region_id']); values=[]
            for sample in sample_ids:
                if not r[sample].isdigit(): raise SystemExit(f'nonnegative integer count required at line {n}, sample {sample}')
                values.append(int(r[sample]))
            base=[r['region_id'],r['chromosome'],s,e]; wo.writerow(base+values); wr.writerow(base+['0-based-half-open'])
    provenance={'route':'externally_processed','source_files':[{'sample_id':r['sample_id'],'source_genome_build':r['source_genome_build'],'source_processing_pipeline':r['source_processing_pipeline'],'source_file_accession':r['source_file_accession'],'source_file_checksum':r['source_file_checksum'],'coordinate_normalization_status':r['coordinate_normalization_status']} for r in meta]}
    open(a.provenance,'w').write(json.dumps(provenance,sort_keys=True,indent=2)+'\n')
if __name__=='__main__': main()
