#!/usr/bin/env python3
import argparse,csv,datetime,hashlib,json
from pathlib import Path

def sha256(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()
def read_tsv(p):
    with open(p,encoding='utf-8-sig',newline='') as h: return list(csv.DictReader(h,delimiter='\t'))
def main():
    p=argparse.ArgumentParser()
    for x in ('samplesheet','reference','differential-status','differential-design','containers','nextflow-version','project-version','profile','design','contrast','run-de','reference-name','reference-scope','fasta-source-url','gtf-source-url','manifest','statuses'): p.add_argument('--'+x,required=True)
    p.add_argument('--versions',nargs='+',required=True); a=p.parse_args()
    rows=read_tsv(a.samplesheet); ref=json.load(open(a.reference,encoding='utf-8')); de=json.load(open(a.differential_status,encoding='utf-8')); design=json.load(open(a.differential_design,encoding='utf-8'))
    containers=read_tsv(a.containers); versions={}
    for f in sorted(a.versions):
        key=Path(f).name
        lines=Path(f).read_text(encoding='utf-8',errors='replace').strip().splitlines()
        versions.setdefault(key,[]).extend(lines if lines else ['UNRESOLVED'])
    input_files=[]
    for row in rows:
        item={'sample_id':row['sample_id'],'library_id':row['library_id'],'route':row['input_route']}
        field='read1' if row['input_route']=='raw' else 'processed_input'; q=Path(row[field])
        item[field]={'name':q.name,'sha256':sha256(q)}
        if row.get('read2'):
            q2=Path(row['read2']); item['read2']={'name':q2.name,'sha256':sha256(q2)}
        input_files.append(item)
    public_meta=[]
    excluded={'read1','read2','processed_input'}
    for row in rows: public_meta.append({k:v for k,v in row.items() if k not in excluded})
    manifest={'project':'regulatory-multiomics-pipeline','project_version':a.project_version,'run_timestamp_utc':datetime.datetime.now(datetime.timezone.utc).replace(microsecond=0).isoformat(),'nextflow_version':a.nextflow_version,'execution_profile':a.profile,'parameters':{'design_formula':a.design,'contrast':a.contrast or 'NOT_PROVIDED','run_de':a.run_de.lower()=='true'},'reference':{'name':a.reference_name,'scope':a.reference_scope,'fasta_source_url':a.fasta_source_url,'gtf_source_url':a.gtf_source_url,'validation':ref},'inputs':input_files,'sample_design':public_meta,'differential_design':design,'analysis_status':de,'tool_versions':versions,'containers':containers}
    Path(a.manifest).write_text(json.dumps(manifest,sort_keys=True,indent=2)+'\n',encoding='utf-8')
    raw_count=sum(r['input_route']=='raw' for r in rows); proc_count=len(rows)-raw_count
    statuses=[('rna_raw_processing','COMPLETE' if raw_count else 'NOT_AVAILABLE_INPUT'),('rna_processed_input','COMPLETE' if proc_count else 'NOT_AVAILABLE_INPUT'),('rna_count_contract','COMPLETE'),('rna_differential_expression',de['status']),('rna_multiqc','COMPLETE'),('rna_provenance','COMPLETE')]
    with open(a.statuses,'w',newline='',encoding='utf-8') as h:
        w=csv.writer(h,delimiter='\t',lineterminator='\n'); w.writerow(['component','status']); w.writerows(statuses)
if __name__=='__main__': main()
