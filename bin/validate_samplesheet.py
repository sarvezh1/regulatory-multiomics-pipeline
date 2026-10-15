#!/usr/bin/env python3
import argparse, csv, json, re
from pathlib import Path

REQUIRED = ['sample_id','library_id','condition','assay','biological_replicate','library_layout','strandedness','read1','read2','processed_input','processed_input_type','source_genome_build','source_annotation','source_processing_pipeline','source_file_accession','source_file_checksum','coordinate_normalization_status','gene_id_mapping_status','genome_build']
SAFE = re.compile(r'^[A-Za-z0-9][A-Za-z0-9_.-]*$')
RNA_NAMES = {'RNA':'RNA','RNA-seq':'RNA','rna':'RNA','rna-seq':'RNA'}
UNKNOWN = {'UNKNOWN','UNRESOLVED'}

def fail(errors):
    raise SystemExit('Samplesheet validation failed:\n- ' + '\n- '.join(errors))

def main():
    p=argparse.ArgumentParser(); p.add_argument('--input',required=True); p.add_argument('--output',required=True); p.add_argument('--report',required=True); a=p.parse_args()
    src=Path(a.input).resolve(); errors=[]
    with src.open(newline='',encoding='utf-8-sig') as h:
        reader=csv.DictReader(h,delimiter='\t')
        if not reader.fieldnames: fail(['empty file or missing header'])
        missing=[x for x in REQUIRED if x not in reader.fieldnames]
        if missing: fail(['missing required columns: '+', '.join(missing)])
        rows=list(reader)
    if not rows: fail(['no data rows'])
    seen_lib=set(); sample_meta={}; normalized=[]; has_raw=False
    for n,row in enumerate(rows,2):
        row={k:(v or '').strip() for k,v in row.items()}
        prefix=f'row {n}'
        if any("'" in v or '\x00' in v for v in row.values()): errors.append(f'{prefix}: apostrophes and NUL characters are not accepted in metadata values')
        for key in ('sample_id','library_id','condition'):
            if not row[key] or not SAFE.fullmatch(row[key]): errors.append(f'{prefix}: unsafe or empty {key}')
        if row['library_id'] in seen_lib: errors.append(f'{prefix}: duplicate library_id {row["library_id"]}')
        seen_lib.add(row['library_id'])
        if row['assay'] not in RNA_NAMES:
            errors.append(f'{prefix}: unsupported assay {row["assay"]!r}; the RNA branch accepts RNA or RNA-seq')
        else: row['assay']='RNA'
        for key in ('biological_replicate','technical_replicate'):
            if key=='technical_replicate' and not row.get(key): row[key]='1'
            if not re.fullmatch(r'[1-9][0-9]*',row.get(key,'')): errors.append(f'{prefix}: {key} must be a positive integer')
        if row['library_layout'] not in {'single','paired','not_applicable'}: errors.append(f'{prefix}: invalid library_layout')
        if row['strandedness'] not in {'unstranded','forward','reverse'}: errors.append(f'{prefix}: strandedness must be explicit')
        raw=bool(row['read1'] or row['read2']); processed=bool(row['processed_input'])
        if raw == processed: errors.append(f'{prefix}: supply exactly one of raw reads or processed_input')
        if raw:
            has_raw=True; row['input_route']='raw'
            if row['library_layout'] not in {'single','paired'}: errors.append(f'{prefix}: raw input requires single or paired layout')
            if not row['read1']: errors.append(f'{prefix}: raw input requires read1')
            if row['library_layout']=='paired' and not row['read2']: errors.append(f'{prefix}: paired input requires read2')
            if row['library_layout']=='single' and row['read2']: errors.append(f'{prefix}: single input forbids read2')
            if row['processed_input_type']: errors.append(f'{prefix}: raw input forbids processed_input_type')
            path_fields=['read1'] + (['read2'] if row['read2'] else [])
        else:
            row['input_route']='processed'; path_fields=['processed_input']
            if row['library_layout']!='not_applicable': errors.append(f'{prefix}: processed input requires not_applicable layout')
            if row['processed_input_type']!='gene_count_tsv': errors.append(f'{prefix}: processed_input_type must be gene_count_tsv')
            for key in ('source_genome_build','source_annotation','source_processing_pipeline','source_file_accession','source_file_checksum','coordinate_normalization_status','gene_id_mapping_status'):
                if not row[key]: errors.append(f'{prefix}: {key} must be recorded or explicitly UNKNOWN/UNRESOLVED')
            if row['source_file_checksum'] not in UNKNOWN|{'NOT_APPLICABLE'} and not re.fullmatch(r'[A-Za-z0-9_-]+:[0-9A-Fa-f]+',row['source_file_checksum']): errors.append(f'{prefix}: source_file_checksum must be algorithm:hex, UNKNOWN, UNRESOLVED, or NOT_APPLICABLE')
            if row['coordinate_normalization_status'] not in {'not_required','verified_primary_assembly','normalized_to_primary','UNKNOWN','UNRESOLVED'}: errors.append(f'{prefix}: invalid coordinate_normalization_status')
            if row['gene_id_mapping_status'] not in {'not_applicable','not_required','mapped_to_gencode50','partially_mapped','UNKNOWN','UNRESOLVED'}: errors.append(f'{prefix}: invalid gene_id_mapping_status')
        for key in path_fields:
            q=Path(row[key]); q=q if q.is_absolute() else src.parent/q
            if not q.is_file(): errors.append(f'{prefix}: {key} does not exist: {row[key]}')
            row[key]=str(q.resolve())
        if row['genome_build']!='GRCh38_primary': errors.append(f'{prefix}: genome_build must be GRCh38_primary')
        stable=tuple(row.get(k,'') for k in ('condition','subject_id','donor_id','biosample_id','biological_replicate','batch','treatment','timepoint','genotype','strandedness','genome_build'))
        if row['sample_id'] in sample_meta and sample_meta[row['sample_id']] != stable:
            errors.append(f'{prefix}: repeated sample_id has contradictory biological metadata')
        sample_meta[row['sample_id']]=stable
        normalized.append(row)
    if errors: fail(errors)
    fields=list(reader.fieldnames)
    if 'technical_replicate' not in fields: fields.append('technical_replicate')
    if 'input_route' not in fields: fields.append('input_route')
    with open(a.output,'w',newline='',encoding='utf-8') as h:
        w=csv.DictWriter(h,fieldnames=fields,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(normalized)
    if has_raw: Path('raw_inputs.present').write_text('raw inputs validated\n',encoding='utf-8')
    report={'status':'COMPLETE','records':len(normalized),'biological_samples':len(sample_meta),'raw_libraries':sum(r['input_route']=='raw' for r in normalized),'processed_libraries':sum(r['input_route']=='processed' for r in normalized)}
    Path(a.report).write_text(json.dumps(report,sort_keys=True,indent=2)+'\n',encoding='utf-8')
if __name__=='__main__': main()
