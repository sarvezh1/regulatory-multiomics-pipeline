#!/usr/bin/env python3
import argparse,csv,hashlib,json
p=argparse.ArgumentParser(); p.add_argument('--matrix',required=True); p.add_argument('--metadata',required=True); p.add_argument('--counts',required=True); p.add_argument('--regions',required=True); p.add_argument('--provenance',required=True); a=p.parse_args()
meta=list(csv.DictReader(open(a.metadata),delimiter='\t')); treatments=[r['sample_id'] for r in meta if r['assay']=='H3K27ac_ChIP']; controls=[r['sample_id'] for r in meta if r['assay']=='ChIP_input']
rows=list(csv.DictReader(open(a.matrix),delimiter='\t')); req=['region_id','chromosome','start','end']+treatments+controls
if not rows or any(x not in (rows[0].keys() if rows else []) for x in req): raise SystemExit('processed H3K27ac matrix lacks coordinates or required treatment/control sample columns')
seen=set()
for n,r in enumerate(rows,2):
    s,e=int(r['start']),int(r['end'])
    if s<0 or e<=s: raise SystemExit(f'invalid coordinates at row {n}')
    key=(r['chromosome'],s,e)
    if key in seen: raise SystemExit(f'duplicate region at row {n}')
    seen.add(key)
    expected='H3K27AC_'+hashlib.sha256(f'{r["chromosome"]}:{s}-{e}'.encode()).hexdigest()[:16]
    if not r['region_id']: r['region_id']=expected
    for x in treatments+controls:
        if not r[x].isdigit(): raise SystemExit(f'{x} count at row {n} is not a nonnegative integer')
with open(a.counts,'w',newline='') as o:
    w=csv.DictWriter(o,fieldnames=req,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows([{k:r[k] for k in req} for r in rows])
with open(a.regions,'w',newline='') as o:
    w=csv.writer(o,delimiter='\t',lineterminator='\n'); w.writerow(['region_id','chromosome','start','end','coordinate_system']); [w.writerow([r['region_id'],r['chromosome'],r['start'],r['end'],'0-based-half-open']) for r in rows]
prov={'external_processed_input':True,'target':'H3K27ac','controls_retained':controls,'sources':[{k:r.get(k,'UNKNOWN') for k in ('sample_id','source_genome_build','source_annotation','source_processing_pipeline','source_file_accession','source_file_checksum','coordinate_normalization_status')} for r in meta]}
open(a.provenance,'w').write(json.dumps(prov,sort_keys=True,indent=2)+'\n')
