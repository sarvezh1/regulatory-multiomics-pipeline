#!/usr/bin/env python3
import argparse, csv, hashlib, json, re
from pathlib import Path

REQUIRED = ['sample_id','library_id','condition','assay','biological_replicate','technical_replicate','library_layout','bisulfite_directionality','bisulfite_protocol','read1','read2','processed_input','processed_input_type','source_genome_build','source_annotation','source_processing_pipeline','source_file_accession','source_file_checksum','coordinate_normalization_status','genome_build']
SAFE = re.compile(r'^[A-Za-z0-9][A-Za-z0-9_.-]*$')
WGBS_NAMES = {'WGBS':'WGBS','wgbs':'WGBS','WGBS-seq':'WGBS','wgbs-seq':'WGBS'}
UNKNOWN = {'UNKNOWN','UNRESOLVED','NOT_APPLICABLE'}

def sha256(path):
    digest=hashlib.sha256()
    with open(path,'rb') as handle:
        for block in iter(lambda:handle.read(1024*1024),b''): digest.update(block)
    return digest.hexdigest()

def fail(errors):
    raise SystemExit('WGBS samplesheet validation failed:\n- '+'\n- '.join(errors))

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--input',required=True); parser.add_argument('--output',required=True); parser.add_argument('--report',required=True)
    args=parser.parse_args(); source=Path(args.input).resolve(); errors=[]
    with source.open(newline='',encoding='utf-8-sig') as handle:
        reader=csv.DictReader(handle,delimiter='\t')
        if not reader.fieldnames: fail(['empty file or missing header'])
        missing=[field for field in REQUIRED if field not in reader.fieldnames]
        if missing: fail(['missing required columns: '+', '.join(missing)])
        fields=list(reader.fieldnames); rows=list(reader)
    if not rows: fail(['no data rows'])
    seen=set(); sample_meta={}; normalized=[]; input_files={}; routes=set(); processed_types=set(); quantitative=True
    for line,raw_row in enumerate(rows,2):
        row={key:(value or '').strip() for key,value in raw_row.items()}; prefix=f'row {line}'
        for key in ('sample_id','library_id','condition'):
            if not row.get(key) or not SAFE.fullmatch(row[key]): errors.append(f'{prefix}: unsafe or empty {key}')
        if row.get('library_id') in seen: errors.append(f'{prefix}: duplicate library_id {row.get("library_id")}')
        seen.add(row.get('library_id'))
        if row.get('assay') not in WGBS_NAMES: errors.append(f'{prefix}: unsupported assay {row.get("assay")!r}; accepted values are WGBS and WGBS-seq')
        else: row['assay']='WGBS'
        if row.get('bisulfite_protocol')!='conventional_wgbs': errors.append(f'{prefix}: bisulfite_protocol must be conventional_wgbs')
        for key in ('biological_replicate','technical_replicate'):
            if not re.fullmatch(r'[1-9][0-9]*',row.get(key,'')): errors.append(f'{prefix}: {key} must be a positive integer')
        raw=bool(row.get('read1') or row.get('read2')); processed=bool(row.get('processed_input'))
        if raw == processed: errors.append(f'{prefix}: supply exactly one of paired raw reads or processed_input')
        if raw:
            route='raw'; paths=('read1','read2')
            if row.get('library_layout')!='paired': errors.append(f'{prefix}: raw WGBS requires library_layout paired')
            if not row.get('read1') or not row.get('read2'): errors.append(f'{prefix}: raw WGBS requires both read1 and read2')
            if row.get('bisulfite_directionality') not in {'directional','non_directional'}: errors.append(f'{prefix}: raw WGBS requires explicit directional or non_directional bisulfite_directionality')
            if row.get('processed_input_type'): errors.append(f'{prefix}: raw input forbids processed_input_type')
        else:
            route='processed'; paths=('processed_input',)
            if row.get('library_layout')!='not_applicable': errors.append(f'{prefix}: processed input requires library_layout not_applicable')
            if row.get('bisulfite_directionality') not in {'directional','non_directional','UNKNOWN','UNRESOLVED'}: errors.append(f'{prefix}: invalid bisulfite_directionality')
            if row.get('processed_input_type') not in {'cpg_count_tsv','cpg_percentage_tsv'}: errors.append(f'{prefix}: processed_input_type must be cpg_count_tsv or cpg_percentage_tsv')
            else: processed_types.add(row.get('processed_input_type'))
            if row.get('processed_input_type')=='cpg_percentage_tsv': quantitative=False
            for key in ('source_genome_build','source_annotation','source_processing_pipeline','source_file_accession','source_file_checksum','coordinate_normalization_status'):
                if not row.get(key): errors.append(f'{prefix}: {key} must be recorded or explicitly UNKNOWN/UNRESOLVED')
            if row.get('source_file_checksum') not in UNKNOWN and not re.fullmatch(r'[A-Za-z0-9_-]+:[0-9A-Fa-f]+',row.get('source_file_checksum','')): errors.append(f'{prefix}: invalid source_file_checksum')
            if row.get('coordinate_normalization_status') not in {'verified_primary_assembly','normalized_to_primary','UNKNOWN','UNRESOLVED'}: errors.append(f'{prefix}: invalid coordinate_normalization_status')
            build=row.get('source_genome_build','').lower()
            if not any(token in build for token in ('grch38','hg38','synthetic')) and row.get('coordinate_normalization_status')!='normalized_to_primary': errors.append(f'{prefix}: incompatible source genome requires normalized_to_primary')
        routes.add(route); row['input_route']=route
        for key in paths:
            value=row.get(key,''); path=Path(value); path=path if path.is_absolute() else source.parent/path
            if not path.is_file(): errors.append(f'{prefix}: {key} does not exist: {value}')
            else: input_files[str(path.resolve())]={'filename':path.name,'size_bytes':path.stat().st_size,'sha256':sha256(path)}
            row[key]=str(path.resolve())
        if row.get('genome_build')!='GRCh38_primary': errors.append(f'{prefix}: genome_build must be GRCh38_primary')
        stable=tuple(row.get(key,'') for key in ('condition','subject_id','donor_id','biosample_id','biological_replicate','batch','treatment','timepoint','genotype','genome_build'))
        if row.get('sample_id') in sample_meta and sample_meta[row['sample_id']] != stable: errors.append(f'{prefix}: repeated sample_id has contradictory biological metadata')
        sample_meta[row.get('sample_id')]=stable; normalized.append(row)
    if len(routes)!=1: errors.append('mixed raw and processed WGBS routes are not supported in one run')
    if len(processed_types)>1: errors.append('mixed count-preserving and percentage-only processed WGBS inputs are not supported in one run')
    if errors: fail(errors)
    if 'input_route' not in fields: fields.append('input_route')
    with open(args.output,'w',newline='',encoding='utf-8') as handle:
        writer=csv.DictWriter(handle,fieldnames=fields,delimiter='\t',lineterminator='\n'); writer.writeheader(); writer.writerows(normalized)
    route=next(iter(routes)); Path(f'{route}_inputs.present').write_text(f'{route} WGBS inputs validated\n')
    report={'status':'COMPLETE','route':route,'records':len(normalized),'biological_samples':len(sample_meta),'quantitative_counts_available':quantitative,'input_files':sorted(input_files.values(),key=lambda item:item['filename'])}
    Path(args.report).write_text(json.dumps(report,sort_keys=True,indent=2)+'\n')

if __name__=='__main__': main()
