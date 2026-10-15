#!/usr/bin/env python3
import argparse,csv,re

def main():
    p=argparse.ArgumentParser()
    for x in ('input','input-type','sample-id','library-id','source-genome-build','source-annotation','source-processing-pipeline','source-file-accession','source-file-checksum','coordinate-normalization-status','gene-id-mapping-status','output'): p.add_argument('--'+x,required=True)
    a=p.parse_args()
    if a.input_type!='gene_count_tsv': raise SystemExit('only gene_count_tsv is supported')
    with open(a.input,encoding='utf-8-sig',newline='') as h: rows=list(csv.DictReader(h,delimiter='\t'))
    required={'gene_id_original','raw_count'}
    if not rows or not required.issubset(rows[0]): raise SystemExit('processed count TSV requires gene_id_original and raw_count')
    seen=set(); fields=['gene_id_original','gene_id_base','gene_name','raw_count','sample_id','library_id','source_genome_build','source_annotation','source_processing_pipeline','source_file_accession','source_file_checksum','coordinate_normalization_status','gene_id_mapping_status']
    with open(a.output,'w',newline='',encoding='utf-8') as h:
        w=csv.DictWriter(h,fields,delimiter='\t',lineterminator='\n'); w.writeheader()
        for i,row in enumerate(rows,2):
            gid=row['gene_id_original'].strip(); value=row['raw_count'].strip()
            if not gid or gid in seen: raise SystemExit(f'empty or duplicate gene ID at row {i}')
            if not re.fullmatch(r'\d+',value): raise SystemExit(f'raw_count must be a non-negative integer at row {i}')
            seen.add(gid); base=row.get('gene_id_base','').strip() or gid.split('.')[0]
            w.writerow({'gene_id_original':gid,'gene_id_base':base,'gene_name':row.get('gene_name','').strip(),'raw_count':value,'sample_id':a.sample_id,'library_id':a.library_id,'source_genome_build':a.source_genome_build,'source_annotation':a.source_annotation,'source_processing_pipeline':a.source_processing_pipeline,'source_file_accession':a.source_file_accession,'source_file_checksum':a.source_file_checksum,'coordinate_normalization_status':a.coordinate_normalization_status,'gene_id_mapping_status':a.gene_id_mapping_status})
if __name__=='__main__': main()
