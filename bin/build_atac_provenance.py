#!/usr/bin/env python3
import argparse, csv, datetime, hashlib, json
from pathlib import Path
def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()
def input_files(samplesheet):
    unique={}
    with open(samplesheet,encoding='utf-8-sig',newline='') as handle:
        for row in csv.DictReader(handle,delimiter='\t'):
            keys=('read1','read2') if row.get('input_route')=='raw' else ('processed_input',)
            for key in keys:
                value=row.get(key,'')
                if not value: continue
                path=Path(value)
                unique[str(path.resolve())]={'filename':path.name,'size_bytes':path.stat().st_size,'sha256':sha(path)}
    return sorted(unique.values(),key=lambda item:item['filename'])
def main():
    p=argparse.ArgumentParser();
    for k in ('samplesheet','validation_report','reference_report','blacklist','tss','status','da_status','versions','containers'): p.add_argument('--'+k.replace('_','-'),required=True)
    p.add_argument('--profile',required=True); p.add_argument('--nextflow-version',required=True); p.add_argument('--design',required=True); p.add_argument('--contrast',default=''); p.add_argument('--effective-genome-size',required=True); p.add_argument('--output',required=True); a=p.parse_args()
    statuses=[]
    for path in (a.status,a.da_status):
        with open(path) as h:
            obj=json.load(h)
            if 'components' in obj: statuses.extend(obj['components'])
            else: statuses.append({'component':'atac_differential','status':obj['status'],'reason':obj.get('reason','')})
    validation=json.load(open(a.validation_report))
    manifest={'project':'regulatory-multiomics-pipeline','version':'0.1.0','run_timestamp_utc':datetime.datetime.now(datetime.timezone.utc).replace(microsecond=0).isoformat(),'nextflow_version':a.nextflow_version,'execution_profile':a.profile,'inputs':input_files(a.samplesheet),'samplesheet':{'sha256':sha(a.samplesheet),'records':validation.get('records'),'biological_samples':validation.get('biological_samples')},'reference':json.load(open(a.reference_report)),'blacklist':{'sha256':sha(a.blacklist),'filename':Path(a.blacklist).name},'tss':{'sha256':sha(a.tss),'filename':Path(a.tss).name,'policy':'deduplicated unique transcript TSS positions'},'parameters':{'bowtie2':'--very-sensitive -X 2000 --no-mixed --no-discordant','mapq_minimum':30,'duplicate_policy':'samtools name-sort, fixmate, coordinate-sort, markdup; duplicate fragments excluded downstream','mitochondrial_policy':'account separately and exclude from nuclear analysis','macs3':'-f BAMPE -q 0.01 --keep-dup all --call-summits','effective_genome_size':int(a.effective_genome_size),'consensus_rule':'overlap-connected merge; book-ended intervals remain separate','counting_rule':'one paired fragment per consensus overlap','design_formula':a.design,'contrast':a.contrast},'statuses':statuses,'tools':Path(a.versions).read_text().splitlines(),'containers':list(csv.DictReader(open(a.containers),delimiter='\t'))}
    Path(a.output).write_text(json.dumps(manifest,sort_keys=True,indent=2)+'\n')
if __name__=='__main__': main()
