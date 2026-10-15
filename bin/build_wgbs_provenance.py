#!/usr/bin/env python3
import argparse, csv, datetime, hashlib, json
from pathlib import Path

def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as handle:
        for block in iter(lambda:handle.read(1024*1024),b''): h.update(block)
    return h.hexdigest()

def inputs(samplesheet):
    found={}
    with open(samplesheet,encoding='utf-8-sig',newline='') as handle:
        for row in csv.DictReader(handle,delimiter='\t'):
            keys=('read1','read2') if row.get('input_route')=='raw' else ('processed_input',)
            for key in keys:
                value=row.get(key,'')
                if not value: continue
                p=Path(value); found[str(p.resolve())]={'filename':p.name,'size_bytes':p.stat().st_size,'sha256':sha(p)}
    return sorted(found.values(),key=lambda x:x['filename'])

def main():
    p=argparse.ArgumentParser()
    for key in ('samplesheet','validation-report','reference-report','reference-checksums','processing-status','dm-status','versions','containers','output'): p.add_argument('--'+key,required=True)
    for key in ('nextflow-version','profile','design','contrast'): p.add_argument('--'+key,required=True)
    p.add_argument('--min-coverage',type=int,required=True); p.add_argument('--ignore',type=int,required=True); p.add_argument('--ignore-r2',type=int,required=True); p.add_argument('--ignore-3prime',type=int,required=True); p.add_argument('--ignore-3prime-r2',type=int,required=True)
    a=p.parse_args(); validation=json.load(open(a.validation_report)); statuses=json.load(open(a.processing_status))['components']; statuses.append(json.load(open(a.dm_status)))
    sample_rows=list(csv.DictReader(open(a.samplesheet,encoding='utf-8'),delimiter='\t'))
    manifest={'project':'regulatory-multiomics-pipeline','version':'0.1.0','run_timestamp_utc':datetime.datetime.now(datetime.timezone.utc).replace(microsecond=0).isoformat(),'nextflow_version':a.nextflow_version,'execution_profile':a.profile,'inputs':inputs(a.samplesheet),'samplesheet':{'sha256':sha(a.samplesheet),'records':validation['records'],'biological_samples':validation['biological_samples'],'libraries':[{'sample_id':r['sample_id'],'library_id':r['library_id'],'bisulfite_protocol':r['bisulfite_protocol'],'bisulfite_directionality':r['bisulfite_directionality'],'conversion_control_contig':r.get('conversion_control_contig','')} for r in sample_rows]},'reference':json.load(open(a.reference_report)),'bismark_reference_checksums':list(csv.DictReader(open(a.reference_checksums),delimiter='\t')) if Path(a.reference_checksums).stat().st_size else [],'parameters':{'trim_galore':'--paired --quality 20 --length 20','bismark_prepare':'prepare --bowtie2','bismark_align':'align --bowtie2 faithful default; --non_directional only when declared','deduplication':'per technical library with bismark dedup --paired','extraction':'paired; no_overlap; comprehensive CpG/CHG/CHH; standard 1-based closed Bismark coverage','mbias_ignore':{'ignore':a.ignore,'ignore_r2':a.ignore_r2,'ignore_3prime':a.ignore_3prime,'ignore_3prime_r2':a.ignore_3prime_r2},'cpg_policy':'reference-aware 1-based source conversion and complementary-strand count summation to 0-based-half-open CpG dinucleotide; non-reference CpG calls are retained in an explicit exclusion table','minimum_cpg_coverage':a.min_coverage,'technical_library_aggregation':'sum methylated and unmethylated counts before fractions','design_formula':a.design,'contrast':a.contrast},'statuses':statuses,'tools':Path(a.versions).read_text().splitlines(),'containers':list(csv.DictReader(open(a.containers),delimiter='\t'))}
    Path(a.output).write_text(json.dumps(manifest,sort_keys=True,indent=2)+'\n')
if __name__=='__main__': main()
