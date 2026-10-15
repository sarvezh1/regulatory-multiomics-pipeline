#!/usr/bin/env python3
import argparse, csv, hashlib, json, re
from pathlib import Path

REQUIRED = ['sample_id','library_id','biosample_id','condition','assay','biological_replicate','technical_replicate','library_layout','read1','read2','processed_input','processed_input_type','control_sample_id','histone_target','source_genome_build','source_annotation','source_processing_pipeline','source_file_accession','source_file_checksum','coordinate_normalization_status','genome_build']
SAFE = re.compile(r'^[A-Za-z0-9][A-Za-z0-9_.-]*$')
UNKNOWN = {'UNKNOWN','UNRESOLVED','NOT_APPLICABLE'}

def digest(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()

def fail(errors): raise SystemExit('H3K27ac samplesheet validation failed:\n- '+'\n- '.join(errors))

def main():
    p=argparse.ArgumentParser(); p.add_argument('--input',required=True); p.add_argument('--output',required=True); p.add_argument('--report',required=True); a=p.parse_args()
    src=Path(a.input).resolve(); errors=[]
    with src.open(newline='',encoding='utf-8-sig') as h:
        r=csv.DictReader(h,delimiter='\t'); fields=list(r.fieldnames or []); rows=list(r)
    missing=[x for x in REQUIRED if x not in fields]
    if missing: fail(['missing required columns: '+', '.join(missing)])
    if not rows: fail(['no data rows'])
    normalized=[]; seen=set(); samples={}; files={}; has_raw=False; has_processed=False
    for n,rawrow in enumerate(rows,2):
        row={k:(v or '').strip() for k,v in rawrow.items()}; pre=f'row {n}'
        for k in ('sample_id','library_id','condition'):
            if not SAFE.fullmatch(row.get(k,'')): errors.append(f'{pre}: unsafe or empty {k}')
        if row['library_id'] in seen: errors.append(f'{pre}: duplicate library_id {row["library_id"]}')
        seen.add(row['library_id'])
        if row['assay'] not in {'H3K27ac_ChIP','ChIP_input'}: errors.append(f'{pre}: assay must be H3K27ac_ChIP or ChIP_input')
        target='H3K27ac' if row['assay']=='H3K27ac_ChIP' else 'CONTROL'
        if row['histone_target']!=target: errors.append(f'{pre}: histone_target must be {target} for {row["assay"]}')
        for k in ('biological_replicate','technical_replicate'):
            if not re.fullmatch(r'[1-9][0-9]*',row.get(k,'')): errors.append(f'{pre}: {k} must be a positive integer')
        israw=bool(row['read1'] or row['read2']); isproc=bool(row['processed_input'])
        if israw==isproc: errors.append(f'{pre}: supply exactly one raw or processed input route')
        paths=[]
        if israw:
            has_raw=True; row['input_route']='raw'; paths=['read1']
            if row['library_layout'] not in {'single','paired'}: errors.append(f'{pre}: raw layout must be single or paired')
            if not row['read1']: errors.append(f'{pre}: raw route requires read1')
            if row['library_layout']=='paired' and not row['read2']: errors.append(f'{pre}: paired layout requires read2')
            if row['library_layout']=='single' and row['read2']: errors.append(f'{pre}: single layout forbids read2')
            if row['read2']: paths.append('read2')
            if row['processed_input_type']: errors.append(f'{pre}: raw input forbids processed_input_type')
        else:
            has_processed=True; row['input_route']='processed'; paths=['processed_input']
            if row['library_layout']!='not_applicable': errors.append(f'{pre}: processed input requires not_applicable layout')
            if row['processed_input_type']!='h3k27ac_region_count_matrix_tsv': errors.append(f'{pre}: processed_input_type must be h3k27ac_region_count_matrix_tsv')
            for k in ('source_genome_build','source_annotation','source_processing_pipeline','source_file_accession','source_file_checksum','coordinate_normalization_status'):
                if not row[k]: errors.append(f'{pre}: {k} must be recorded or explicitly UNKNOWN/UNRESOLVED')
            if row['source_file_checksum'] not in UNKNOWN and not re.fullmatch(r'[A-Za-z0-9_-]+:[0-9A-Fa-f]+',row['source_file_checksum']): errors.append(f'{pre}: invalid source_file_checksum')
            source=row['source_genome_build'].lower()
            if not any(x in source for x in ('grch38','hg38','synthetic')) and row['coordinate_normalization_status']!='normalized_to_primary': errors.append(f'{pre}: incompatible source genome without normalized_to_primary')
        for k in paths:
            q=Path(row[k]); q=q if q.is_absolute() else src.parent/q
            if not q.is_file(): errors.append(f'{pre}: {k} does not exist: {row[k]}')
            else: files[str(q.resolve())]={'filename':q.name,'size_bytes':q.stat().st_size,'sha256':digest(q)}
            row[k]=str(q.resolve())
        if row['genome_build']!='GRCh38_primary': errors.append(f'{pre}: genome_build must be GRCh38_primary')
        stable=tuple(row.get(k,'') for k in ('assay','condition','biosample_id','biological_replicate','batch','treatment','timepoint','genotype','library_layout'))
        if row['sample_id'] in samples and samples[row['sample_id']]!=stable: errors.append(f'{pre}: repeated sample_id has contradictory biological metadata')
        samples[row['sample_id']]=stable; normalized.append(row)
    if has_raw and has_processed: errors.append('mixed raw and processed H3K27ac routes are not supported in one run')
    by_sample={r['sample_id']:r for r in normalized}
    for r in normalized:
        if r['assay']=='H3K27ac_ChIP':
            cid=r['control_sample_id']
            if not cid or cid not in by_sample: errors.append(f'{r["sample_id"]}: control_sample_id does not reference a control record')
            elif by_sample[cid]['assay']!='ChIP_input': errors.append(f'{r["sample_id"]}: control_sample_id must reference ChIP_input')
            else:
                c=by_sample[cid]
                for k in ('condition','biosample_id','library_layout'):
                    if r.get(k)!=c.get(k): errors.append(f'{r["sample_id"]}: control {cid} has incompatible {k}')
        elif r['control_sample_id']: errors.append(f'{r["sample_id"]}: control rows must not have control_sample_id')
    if has_processed:
        mats={r['processed_input'] for r in normalized}
        if len(mats)!=1: errors.append('processed rows must reference one common region count matrix')
    if errors: fail(errors)
    if 'input_route' not in fields: fields.append('input_route')
    with open(a.output,'w',newline='',encoding='utf-8') as h:
        w=csv.DictWriter(h,fieldnames=fields,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(normalized)
    Path('raw_inputs.present' if has_raw else 'processed_inputs.present').write_text('H3K27ac inputs validated\n')
    report={'status':'COMPLETE','records':len(normalized),'biological_samples':len(samples),'treatment_samples':len({r['sample_id'] for r in normalized if r['assay']=='H3K27ac_ChIP'}),'control_samples':len({r['sample_id'] for r in normalized if r['assay']=='ChIP_input'}),'input_files':sorted(files.values(),key=lambda x:x['filename'])}
    Path(a.report).write_text(json.dumps(report,sort_keys=True,indent=2)+'\n')
if __name__=='__main__': main()
