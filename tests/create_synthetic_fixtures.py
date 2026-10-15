#!/usr/bin/env python3
import csv,gzip,random,shutil
from pathlib import Path

COLUMNS=['sample_id','library_id','subject_id','donor_id','biosample_id','condition','assay','biological_replicate','technical_replicate','batch','sex','treatment','timepoint','genotype','library_layout','strandedness','read1','read2','processed_input','processed_input_type','source_genome_build','source_annotation','source_processing_pipeline','source_file_accession','source_file_checksum','coordinate_normalization_status','gene_id_mapping_status','control_sample_id','histone_target','genome_build']
def write_tsv(path,rows,fields=COLUMNS):
    with open(path,'w',newline='',encoding='utf-8') as h:
        w=csv.DictWriter(h,fields,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(rows)
def base(**kw):
    r={x:'' for x in COLUMNS}; r.update({'subject_id':'NA','donor_id':'NA','biosample_id':'NA','assay':'RNA','technical_replicate':'1','batch':'batch1','sex':'not_reported','treatment':'none','timepoint':'0h','genotype':'wild_type','strandedness':'unstranded','control_sample_id':'','histone_target':'','genome_build':'GRCh38_primary'}); r.update(kw); return r
def write_minimal_fastq(path,name):
    path.parent.mkdir(parents=True,exist_ok=True)
    record=f'@{name}\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n'.encode()
    with path.open('wb') as raw:
        with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0) as compressed:
            compressed.write(record)
def write_contract_raw_fixtures():
    fixtures=Path(__file__).parent/'fixtures'
    for relative,name in (
        ('atac/raw/atac_lib_R1.fastq.gz','atac_lib/1'),
        ('atac/raw/atac_lib_R2.fastq.gz','atac_lib/2'),
        ('chip/raw/se_chip.fastq.gz','se_chip'),
        ('chip/raw/se_input.fastq.gz','se_input'),
        ('chip/raw/pe_chip_R1.fastq.gz','pe_chip/1'),
        ('chip/raw/pe_chip_R2.fastq.gz','pe_chip/2'),
        ('chip/raw/pe_input_R1.fastq.gz','pe_input/1'),
        ('chip/raw/pe_input_R2.fastq.gz','pe_input/2'),
    ):
        write_minimal_fastq(fixtures/relative,name)
def main():
    root=Path(__file__).parent/'fixtures'/'generated'
    if root.exists(): shutil.rmtree(root)
    (root/'raw').mkdir(parents=True); (root/'processed').mkdir()
    write_contract_raw_fixtures()
    random.seed(417); seq=''.join(random.choice('ACGT') for _ in range(3000))
    (root/'reference.fa').write_text('>chr22\n'+'\n'.join(seq[i:i+60] for i in range(0,len(seq),60))+'\n')
    (root/'reference.gtf').write_text('## synthetic execution fixture\nchr22\tfixture\tgene\t101\t500\t.\t+\t.\tgene_id "ENSGTEST000001.1"; gene_name "UPTEST";\nchr22\tfixture\texon\t101\t500\t.\t+\t.\tgene_id "ENSGTEST000001.1"; gene_name "UPTEST"; transcript_id "ENSTTEST000001.1";\nchr22\tfixture\tgene\t1501\t1900\t.\t-\t.\tgene_id "ENSGTEST000002.1"; gene_name "DOWNTEST";\nchr22\tfixture\texon\t1501\t1900\t.\t-\t.\tgene_id "ENSGTEST000002.1"; gene_name "DOWNTEST"; transcript_id "ENSTTEST000002.1";\n')
    fq=root/'raw'/'raw_lib.fastq.gz'
    with gzip.open(fq,'wt',newline='\n') as h:
        for i,start in enumerate(list(range(110,350,20))+list(range(1510,1750,20)),1):
            read=seq[start:start+50]; h.write(f'@read{i}\n{read}\n+\n'+('I'*50)+'\n')
    raw=base(sample_id='raw_sample',library_id='raw_lib',condition='control',biological_replicate='1',library_layout='single',read1='raw/raw_lib.fastq.gz')
    write_tsv(root/'raw_samplesheet.tsv',[raw])
    genes=[('ENSGTEST%06d.1'%i,'GENE%d'%i) for i in range(1,101)]
    rows=[]
    for group in ('control','treatment'):
        for rep in range(1,4):
            sample=f'{group}_{rep}'; lib=f'{sample}_lib'; path=root/'processed'/f'{lib}.tsv'
            with open(path,'w',newline='',encoding='utf-8') as h:
                w=csv.writer(h,delimiter='\t',lineterminator='\n'); w.writerow(['gene_id_original','gene_name','raw_count'])
                for i,(gid,name) in enumerate(genes,1):
                    if i==100: value=0
                    elif i>=95: value=(i+rep)%3
                    else:
                        baseline=40+i*7
                        jitter=((i*rep*13)%31)-15
                        value=max(1,baseline+jitter)
                        if 2<=i<=10 and group=='treatment': value*=8
                        if 11<=i<=19 and group=='control': value*=8
                    w.writerow([gid,name,value])
            rows.append(base(sample_id=sample,library_id=lib,condition=group,biological_replicate=str(rep),library_layout='not_applicable',processed_input=f'processed/{lib}.tsv',processed_input_type='gene_count_tsv',source_genome_build='synthetic',source_annotation='synthetic_v1',source_processing_pipeline='fixture_counts_v1',source_file_accession='NOT_APPLICABLE',source_file_checksum='UNRESOLVED',coordinate_normalization_status='not_required',gene_id_mapping_status='not_required'))
    write_tsv(root/'processed_samplesheet.tsv',rows)
    write_tsv(root/'one_vs_one.tsv',[rows[0],rows[3]])
    tech_path=root/'processed'/'control_1_lib_tech2.tsv'
    with open(tech_path,'w',newline='',encoding='utf-8') as h:
        w=csv.writer(h,delimiter='\t',lineterminator='\n'); w.writerow(['gene_id_original','gene_name','raw_count'])
        for gid,name in genes: w.writerow([gid,name,5])
    tech=base(**{**rows[0],'library_id':'control_1_lib_tech2','technical_replicate':'2','processed_input':'processed/control_1_lib_tech2.tsv'})
    write_tsv(root/'technical_replicates.tsv',rows+[tech])
    conflict=base(sample_id='bad',library_id='bad_lib',condition='control',biological_replicate='1',library_layout='single',read1='raw/raw_lib.fastq.gz',processed_input='processed/control_1_lib.tsv',processed_input_type='gene_count_tsv')
    write_tsv(root/'raw_processed_conflict.tsv',[conflict])
    badpair=base(sample_id='badpair',library_id='badpair_lib',condition='control',biological_replicate='1',library_layout='paired',read1='raw/raw_lib.fastq.gz')
    write_tsv(root/'malformed_paired.tsv',[badpair])
    print(root)
if __name__=='__main__': main()
