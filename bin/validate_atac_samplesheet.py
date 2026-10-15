#!/usr/bin/env python3
import argparse, csv, hashlib, json, re
from pathlib import Path

REQUIRED = ['sample_id','library_id','condition','assay','biological_replicate','technical_replicate','library_layout','read1','read2','processed_input','processed_input_type','source_genome_build','source_processing_pipeline','source_file_accession','source_file_checksum','coordinate_normalization_status','genome_build']
SAFE = re.compile(r'^[A-Za-z0-9][A-Za-z0-9_.-]*$')
ATAC_NAMES = {'ATAC':'ATAC','ATAC-seq':'ATAC','atac':'ATAC','atac-seq':'ATAC'}
UNKNOWN = {'UNKNOWN','UNRESOLVED','NOT_APPLICABLE'}

def die(errors):
    raise SystemExit('ATAC samplesheet validation failed:\n- ' + '\n- '.join(errors))

def sha256(path):
    digest=hashlib.sha256()
    with open(path,'rb') as handle:
        for block in iter(lambda:handle.read(1024*1024),b''): digest.update(block)
    return digest.hexdigest()

def main():
    p=argparse.ArgumentParser(); p.add_argument('--input',required=True); p.add_argument('--output',required=True); p.add_argument('--report',required=True); a=p.parse_args()
    src=Path(a.input).resolve(); errors=[]
    with src.open(newline='',encoding='utf-8-sig') as h:
        reader=csv.DictReader(h,delimiter='\t')
        if not reader.fieldnames: die(['empty file or missing header'])
        missing=[x for x in REQUIRED if x not in reader.fieldnames]
        if missing: die(['missing required columns: '+', '.join(missing)])
        rows=list(reader); fields=list(reader.fieldnames)
    if not rows: die(['no data rows'])
    seen=set(); sample_meta={}; normalized=[]; has_raw=False; has_processed=False; input_files={}
    for n,row in enumerate(rows,2):
        row={k:(v or '').strip() for k,v in row.items()}; prefix=f'row {n}'
        for key in ('sample_id','library_id','condition'):
            if not row.get(key) or not SAFE.fullmatch(row[key]): errors.append(f'{prefix}: unsafe or empty {key}')
        if row.get('library_id') in seen: errors.append(f'{prefix}: duplicate library_id {row.get("library_id")}')
        seen.add(row.get('library_id'))
        if row.get('assay') not in ATAC_NAMES: errors.append(f'{prefix}: unsupported assay {row.get("assay")!r}; accepted values are ATAC and ATAC-seq')
        else: row['assay']='ATAC'
        for key in ('biological_replicate','technical_replicate'):
            if not re.fullmatch(r'[1-9][0-9]*',row.get(key,'')): errors.append(f'{prefix}: {key} must be a positive integer')
        raw=bool(row.get('read1') or row.get('read2')); processed=bool(row.get('processed_input'))
        if raw == processed: errors.append(f'{prefix}: supply exactly one of paired raw reads or processed_input')
        if raw:
            has_raw=True; row['input_route']='raw'
            if row.get('library_layout')!='paired': errors.append(f'{prefix}: raw ATAC requires library_layout paired')
            if not row.get('read1') or not row.get('read2'): errors.append(f'{prefix}: raw ATAC requires both read1 and read2')
            if row.get('processed_input_type'): errors.append(f'{prefix}: raw input forbids processed_input_type')
            paths=['read1','read2']
        else:
            has_processed=True; row['input_route']='processed'; paths=['processed_input']
            if row.get('library_layout')!='not_applicable': errors.append(f'{prefix}: processed input requires library_layout not_applicable')
            if row.get('processed_input_type')!='region_count_matrix_tsv': errors.append(f'{prefix}: processed_input_type must be region_count_matrix_tsv')
            for key in ('source_genome_build','source_processing_pipeline','source_file_accession','source_file_checksum','coordinate_normalization_status'):
                if not row.get(key): errors.append(f'{prefix}: {key} must be recorded or explicitly UNKNOWN/UNRESOLVED')
            if row.get('source_file_checksum') not in UNKNOWN and not re.fullmatch(r'[A-Za-z0-9_-]+:[0-9A-Fa-f]+',row.get('source_file_checksum','')): errors.append(f'{prefix}: invalid source_file_checksum')
            if row.get('coordinate_normalization_status') not in {'verified_primary_assembly','normalized_to_primary','UNKNOWN','UNRESOLVED'}: errors.append(f'{prefix}: invalid coordinate_normalization_status')
            source=row.get('source_genome_build','').lower()
            if not any(x in source for x in ('grch38','hg38','synthetic')) and row.get('coordinate_normalization_status')!='normalized_to_primary': errors.append(f'{prefix}: incompatible source genome requires explicit normalized_to_primary status')
        for key in paths:
            value=row.get(key,''); q=Path(value); q=q if q.is_absolute() else src.parent/q
            if not q.is_file(): errors.append(f'{prefix}: {key} does not exist: {value}')
            else: input_files[str(q.resolve())]={'filename':q.name,'size_bytes':q.stat().st_size,'sha256':sha256(q)}
            row[key]=str(q.resolve())
        if row.get('genome_build')!='GRCh38_primary': errors.append(f'{prefix}: genome_build must be GRCh38_primary')
        stable=tuple(row.get(k,'') for k in ('condition','subject_id','donor_id','biosample_id','biological_replicate','batch','treatment','timepoint','genotype','genome_build'))
        if row.get('sample_id') in sample_meta and sample_meta[row['sample_id']] != stable: errors.append(f'{prefix}: repeated sample_id has contradictory biological metadata')
        sample_meta[row.get('sample_id')]=stable; normalized.append(row)
    if has_raw and has_processed: errors.append('mixed raw and processed ATAC routes are not supported in one ATAC run')
    if has_processed:
        matrices={r['processed_input'] for r in normalized}
        if len(matrices)!=1: errors.append('processed ATAC rows must reference one common region count matrix')
    if errors: die(errors)
    if 'input_route' not in fields: fields.append('input_route')
    with open(a.output,'w',newline='',encoding='utf-8') as h:
        w=csv.DictWriter(h,fieldnames=fields,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(normalized)
    if has_raw: Path('raw_inputs.present').write_text('raw ATAC inputs validated\n')
    if has_processed: Path('processed_inputs.present').write_text('processed ATAC inputs validated\n')
    report={'status':'COMPLETE','records':len(normalized),'biological_samples':len(sample_meta),'raw_libraries':sum(r['input_route']=='raw' for r in normalized),'processed_rows':sum(r['input_route']=='processed' for r in normalized),'input_files':sorted(input_files.values(),key=lambda x:x['filename'])}
    Path(a.report).write_text(json.dumps(report,sort_keys=True,indent=2)+'\n')
if __name__=='__main__': main()
