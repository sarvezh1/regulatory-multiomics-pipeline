#!/usr/bin/env python3
import argparse,csv,re
from pathlib import Path

def attrs(s):
    return {m.group(1):m.group(2) for m in re.finditer(r'(\S+) "([^"]*)";',s)}
def main():
    p=argparse.ArgumentParser()
    for x in ('counts','gtf','sample-id','library-id','source-genome-build','source-annotation','source-processing-pipeline','source-file-accession','source-file-checksum','coordinate-normalization-status','gene-id-mapping-status','output'): p.add_argument('--'+x,required=True)
    a=p.parse_args(); ann={}
    with open(a.gtf,encoding='utf-8') as h:
        for line in h:
            if line.startswith('#'): continue
            f=line.rstrip('\n').split('\t')
            if len(f)==9 and f[2]=='gene':
                z=attrs(f[8]); gid=z.get('gene_id');
                if gid: ann[gid]=z.get('gene_name','')
    with open(a.counts,encoding='utf-8') as h:
        lines=(x for x in h if not x.startswith('#')); r=csv.DictReader(lines,delimiter='\t'); rows=list(r)
    count_col=list(rows[0])[-1] if rows else None
    fields=['gene_id_original','gene_id_base','gene_name','raw_count','sample_id','library_id','source_genome_build','source_annotation','source_processing_pipeline','source_file_accession','source_file_checksum','coordinate_normalization_status','gene_id_mapping_status']
    with open(a.output,'w',newline='',encoding='utf-8') as h:
        w=csv.DictWriter(h,fields,delimiter='\t',lineterminator='\n'); w.writeheader()
        for row in rows:
            gid=row['Geneid']; count=row[count_col]
            if not re.fullmatch(r'\d+',count): raise SystemExit(f'non-integer featureCounts value for {gid}')
            w.writerow({'gene_id_original':gid,'gene_id_base':gid.split('.')[0],'gene_name':ann.get(gid,''),'raw_count':count,'sample_id':a.sample_id,'library_id':a.library_id,'source_genome_build':a.source_genome_build,'source_annotation':a.source_annotation,'source_processing_pipeline':a.source_processing_pipeline,'source_file_accession':a.source_file_accession,'source_file_checksum':a.source_file_checksum,'coordinate_normalization_status':a.coordinate_normalization_status,'gene_id_mapping_status':a.gene_id_mapping_status})
if __name__=='__main__': main()
