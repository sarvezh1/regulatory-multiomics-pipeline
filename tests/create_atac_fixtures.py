#!/usr/bin/env python3
import csv, gzip, random, shutil
from pathlib import Path

COLUMNS=['sample_id','library_id','subject_id','donor_id','biosample_id','condition','assay','biological_replicate','technical_replicate','batch','sex','treatment','timepoint','genotype','library_layout','strandedness','read1','read2','processed_input','processed_input_type','source_genome_build','source_annotation','source_processing_pipeline','source_file_accession','source_file_checksum','coordinate_normalization_status','gene_id_mapping_status','control_sample_id','histone_target','genome_build']
def rc(s): return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
def write_tsv(path,rows):
    with open(path,'w',newline='') as h: w=csv.DictWriter(h,COLUMNS,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(rows)
def base(**kw):
    r={x:'' for x in COLUMNS}; r.update(subject_id='NA',donor_id='NA',biosample_id='NA',assay='ATAC',technical_replicate='1',batch='batch1',sex='not_reported',treatment='none',timepoint='0h',genotype='wild_type',strandedness='unstranded',genome_build='GRCh38_primary'); r.update(kw); return r
def main():
    root=Path(__file__).parent/'fixtures'/'atac';
    if root.exists(): shutil.rmtree(root)
    (root/'raw').mkdir(parents=True)
    rng=random.Random(8221); chr1=''.join(rng.choice('ACGT') for _ in range(100000)); chrm=''.join(rng.choice('ACGT') for _ in range(10000))
    with open(root/'reference.fa','w') as h:
        for name,seq in [('chr1',chr1),('chrM',chrm)]: h.write(f'>{name}\n'+'\n'.join(seq[i:i+60] for i in range(0,len(seq),60))+'\n')
    (root/'reference.gtf').write_text('## synthetic ATAC execution fixture\nchr1\tfixture\tgene\t20001\t22000\t.\t+\t.\tgene_id "GENE1.1"; gene_name "OPEN1";\nchr1\tfixture\ttranscript\t20001\t22000\t.\t+\t.\tgene_id "GENE1.1"; transcript_id "TX1.1";\nchr1\tfixture\tgene\t68001\t70000\t.\t-\t.\tgene_id "GENE2.1"; gene_name "OPEN2";\nchr1\tfixture\ttranscript\t68001\t70000\t.\t-\t.\tgene_id "GENE2.1"; transcript_id "TX2.1";\n')
    (root/'blacklist.bed').write_text('chr1\t39800\t40600\tsynthetic_blacklist\n')
    pairs=[]
    def add(chrom,seq,start,length,label): pairs.append((chrom,start,length,label,seq[start:start+50],rc(seq[start+length-50:start+length])))
    for i in range(350): add('chr1',chr1,19800+(i*17)%500,80+(i%70),f'peak1_{i}')
    for i in range(350): add('chr1',chr1,69600+(i*19)%500,90+(i%100),f'peak2_{i}')
    for i in range(100): add('chr1',chr1,1000+(i*701)%95000,180+(i%180),f'background_{i}')
    for i in range(60): add('chr1',chr1,40000+(i%80),100,f'blacklist_{i}')
    for i in range(50): add('chrM',chrm,500+(i*97)%8000,100,f'mt_{i}')
    for i in range(40): add('chr1',chr1,20050,100,f'duplicate_{i}')
    with gzip.open(root/'raw'/'atac_lib_R1.fastq.gz','wt',newline='\n') as r1, gzip.open(root/'raw'/'atac_lib_R2.fastq.gz','wt',newline='\n') as r2:
        for i,(_,_,_,label,s1,s2) in enumerate(pairs,1):
            q='I'*50; r1.write(f'@{label}_{i}/1\n{s1}\n+\n{q}\n'); r2.write(f'@{label}_{i}/2\n{s2}\n+\n{q}\n')
    raw=base(sample_id='atac_sample',library_id='atac_lib',condition='control',biological_replicate='1',library_layout='paired',read1='raw/atac_lib_R1.fastq.gz',read2='raw/atac_lib_R2.fastq.gz')
    write_tsv(root/'raw_samplesheet.tsv',[raw])
    samples=[]
    for group in ('control','treatment'):
        for rep in range(1,4): samples.append(base(sample_id=f'{group}_{rep}',library_id=f'{group}_{rep}_lib',condition=group,biological_replicate=str(rep),library_layout='not_applicable',processed_input='processed_atac_counts.tsv',processed_input_type='region_count_matrix_tsv',source_genome_build='synthetic_GRCh38_like',source_processing_pipeline='fixture_v1',source_file_accession='NOT_APPLICABLE',source_file_checksum='UNRESOLVED',coordinate_normalization_status='normalized_to_primary'))
    write_tsv(root/'processed_samplesheet.tsv',samples)
    with open(root/'processed_atac_counts.tsv','w',newline='') as h:
        w=csv.writer(h,delimiter='\t',lineterminator='\n'); ids=[r['sample_id'] for r in samples]; w.writerow(['region_id','chromosome','start','end']+ids)
        for i in range(1,101):
            c=[50+(i*j)%13 for j in range(1,7)]
            if i<=15: c[3:]=[x*8 for x in c[3:]]
            elif i<=30: c[:3]=[x*8 for x in c[:3]]
            elif i>=95: c=[0]*6
            w.writerow([f'ATAC_fixture_{i:04d}','chr22',(i-1)*1000,(i-1)*1000+500]+c)
    write_tsv(root/'one_vs_one.tsv',[samples[0],samples[3]])
    rank=[dict(r, batch=r['condition']) for r in samples]; write_tsv(root/'rank_deficient.tsv',rank)
    technical=samples+[dict(samples[0],library_id='control_1_tech2',technical_replicate='2')]; write_tsv(root/'processed_technical_duplicate.tsv',technical)
    write_tsv(root/'missing_r2.tsv',[base(sample_id='bad',library_id='bad_lib',condition='control',biological_replicate='1',library_layout='paired',read1='raw/atac_lib_R1.fastq.gz')])
    conflict=base(sample_id='bad2',library_id='bad2_lib',condition='control',biological_replicate='1',library_layout='paired',read1='raw/atac_lib_R1.fastq.gz',read2='raw/atac_lib_R2.fastq.gz',processed_input='processed_atac_counts.tsv',processed_input_type='region_count_matrix_tsv'); write_tsv(root/'raw_processed_conflict.tsv',[conflict])
    invalid=base(**{**raw,'assay':'RNA'}); write_tsv(root/'invalid_assay.tsv',[invalid])
    incompatible=base(**{**samples[0],'sample_id':'incompatible','library_id':'incompatible_lib','source_genome_build':'hg19','coordinate_normalization_status':'UNRESOLVED'}); write_tsv(root/'incompatible_genome.tsv',[incompatible])
    (root/'malformed_blacklist.bed').write_text('chr1\t500\t100\tbad\n')
    print(f'ATAC fixture pairs={len(pairs)} path={root}')
if __name__=='__main__': main()
