#!/usr/bin/env python3
import argparse,csv
from collections import defaultdict

def read_tsv(p):
    with open(p,encoding='utf-8-sig',newline='') as h: return list(csv.DictReader(h,delimiter='\t'))
def write_tsv(p,fields,rows):
    with open(p,'w',newline='',encoding='utf-8') as h:
        w=csv.DictWriter(h,fields,delimiter='\t',lineterminator='\n',extrasaction='ignore'); w.writeheader(); w.writerows(rows)
def main():
    p=argparse.ArgumentParser(); p.add_argument('--samplesheet',required=True); p.add_argument('--contracts',nargs='+',required=True); p.add_argument('--matrix',required=True); p.add_argument('--metadata',required=True); p.add_argument('--annotations',required=True); p.add_argument('--technical-manifest',required=True); a=p.parse_args()
    sheet=read_tsv(a.samplesheet); libmeta={r['library_id']:r for r in sheet}; counts=defaultdict(lambda:defaultdict(int)); annotations={}; provenance={}; manifest=[]
    for path in sorted(a.contracts):
        rows=read_tsv(path)
        if not rows: raise SystemExit(f'empty contract: {path}')
        lib=rows[0]['library_id']; sample=rows[0]['sample_id']
        if lib not in libmeta or libmeta[lib]['sample_id']!=sample: raise SystemExit(f'contract metadata mismatch: {path}')
        manifest.append({'sample_id':sample,'library_id':lib,'technical_replicate':libmeta[lib].get('technical_replicate','1'),'combination_rule':'integer_sum_before_DE'})
        for r in rows:
            gid=r['gene_id_original']; key=(r['gene_id_base'],r['gene_name'])
            if gid in annotations and annotations[gid]!=key: raise SystemExit(f'conflicting annotation for {gid}')
            annotations[gid]=key; counts[gid][sample]+=int(r['raw_count'])
            provenance[(gid,lib)]=r
    samples=sorted({r['sample_id'] for r in sheet}); genes=sorted(annotations)
    matrix_rows=[]
    for gid in genes:
        base,name=annotations[gid]; row={'gene_id_original':gid,'gene_id_base':base,'gene_name':name}; row.update({s:counts[gid].get(s,0) for s in samples}); matrix_rows.append(row)
    write_tsv(a.matrix,['gene_id_original','gene_id_base','gene_name']+samples,matrix_rows)
    meta_rows=[]
    for sample in samples:
        group=[r for r in sheet if r['sample_id']==sample]; r=group[0]
        out={k:v for k,v in r.items() if k not in {'read1','read2','processed_input','library_id','technical_replicate','input_route'}}
        out['sample_id']=sample; out['library_count']=str(len(group)); meta_rows.append(out)
    meta_fields=['sample_id']+[x for x in meta_rows[0] if x!='sample_id']
    write_tsv(a.metadata,meta_fields,meta_rows)
    write_tsv(a.annotations,['gene_id_original','gene_id_base','gene_name'],[{'gene_id_original':g,'gene_id_base':annotations[g][0],'gene_name':annotations[g][1]} for g in genes])
    write_tsv(a.technical_manifest,['sample_id','library_id','technical_replicate','combination_rule'],sorted(manifest,key=lambda x:(x['sample_id'],x['library_id'])))
if __name__=='__main__': main()
