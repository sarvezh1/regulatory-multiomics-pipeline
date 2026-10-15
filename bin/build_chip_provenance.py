#!/usr/bin/env python3
import argparse,csv,datetime,hashlib,json
from pathlib import Path
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
p=argparse.ArgumentParser()
for k in ('samplesheet','validation','reference','blacklist','processing_status','diff_status','versions','containers'):p.add_argument('--'+k.replace('_','-'),required=True)
for k in ('nextflow_version','profile','design','contrast','effective_genome_size','mapq','se_validation_extsize'):p.add_argument('--'+k.replace('_','-'),required=True)
p.add_argument('--output',required=True);a=p.parse_args(); statuses=[]
for f in (a.processing_status,a.diff_status):
    d=json.load(open(f))
    if 'components' in d: statuses.extend(d['components'])
    else: statuses.append({'component':'h3k27ac_differential','status':d['status'],'reason':d.get('reason','')})
rows=list(csv.DictReader(open(a.samplesheet),delimiter='\t')); inp=[]
for r in rows:
    for k in ('read1','read2','processed_input'):
        if r.get(k):
            q=Path(r[k]); inp.append({'filename':q.name,'size_bytes':q.stat().st_size,'sha256':sha(q)})
manifest={'project':'regulatory-multiomics-pipeline','version':'0.1.0','author':'Sarvesh Galgale','run_timestamp_utc':datetime.datetime.now(datetime.timezone.utc).replace(microsecond=0).isoformat(),'nextflow_version':a.nextflow_version,'execution_profile':a.profile,'target':'H3K27ac','inputs':sorted(inp,key=lambda x:x['filename']),'control_mapping':sorted([{'treatment_sample_id':r['sample_id'],'control_sample_id':r['control_sample_id']} for r in rows if r['assay']=='H3K27ac_ChIP'],key=lambda x:x['treatment_sample_id']),'reference':json.load(open(a.reference)),'blacklist':{'filename':Path(a.blacklist).name,'sha256':sha(a.blacklist)},'parameters':{'bowtie2_single':'--very-sensitive','bowtie2_paired':'--very-sensitive --no-mixed --no-discordant -X 1000','mapq_minimum':int(a.mapq),'duplicate_policy':'per-library samtools markdup; duplicates excluded downstream','mitochondrial_policy':'account separately and exclude chrM from nuclear analysis','macs3_single':'-f BAM -q 0.01 --keep-dup all --call-summits -B --SPMR; default model unless the recorded compact-validation extension override is nonzero','macs3_paired':'-f BAMPE -q 0.01 --keep-dup all --call-summits -B --SPMR; observed fragments','single_end_compact_validation_extension':int(a.se_validation_extsize),'effective_genome_size':int(a.effective_genome_size),'signal':'MACS3 SPMR treatment pileup and control lambda; bdgcmp FE','FRiP_denominator':'usable nonduplicate nuclear reads/fragments after filtering','consensus':'strict overlap-connected merge; book-ended intervals remain separate','counting':'single-end reads once; paired-end fragments once','design_formula':a.design,'contrast':a.contrast},'statuses':statuses,'tools':Path(a.versions).read_text().splitlines(),'containers':list(csv.DictReader(open(a.containers),delimiter='\t'))}
Path(a.output).write_text(json.dumps(manifest,sort_keys=True,indent=2)+'\n')
