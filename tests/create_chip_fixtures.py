#!/usr/bin/env python3
import csv,gzip,hashlib,random
from pathlib import Path
root=Path(__file__).parent/'fixtures'/'chip'; raw=root/'raw'; proc=root/'processed'; raw.mkdir(parents=True,exist_ok=True);proc.mkdir(parents=True,exist_ok=True)
def wslpath(p):
    q=str(p.resolve()).replace('\\','/')
    return '/mnt/'+q[0].lower()+q[2:] if len(q)>2 and q[1]==':' else q
rng=random.Random(427); seq=''.join(rng.choice('ACGT') for _ in range(2000000)); mt=''.join(rng.choice('ACGT') for _ in range(5000))
(root/'reference.fa').write_text('>chr1\n'+seq+'\n>chrM\n'+mt+'\n');(root/'reference.gtf').write_text('chr1\tfixture\tgene\t1\t2000000\t.\t+\t.\tgene_id "GENE1.1"; gene_name "GENE1";\n');(root/'blacklist.bed').write_text('chr1\t5000\t5500\n')
def rc(s):return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
def write_se(path,kind):
    reads=[]
    if kind=='chip':
        i=0
        for center in range(20000,1820000,15000):
            for j in range(60):
                reads.append((f'{kind}_{i}',seq[center-130+j:center-80+j]));i+=1
                reads.append((f'{kind}_{i}',rc(seq[center+50+j:center+100+j])));i+=1
        for center in (200000,1400000):
            for j in range(250):
                start=center-400+(j*13)%800; s=seq[start:start+50]; reads.append((f'{kind}_{i}',rc(s) if j%2 else s));i+=1
        for j in range(50): reads.append((f'{kind}_{i+j}',seq[5100+j:5150+j]))
    else:
        for i in range(2000):
            start=(i*947)%1990000; reads.append((f'{kind}_{i}',seq[start:start+50]))
    for j in range(20): reads.append((f'{kind}_mt_{j}',mt[100+j:150+j]))
    reads.extend((f'{n}_duplicate_{i}',s) for i,(n,s) in enumerate(reads[:20]))
    with gzip.open(path,'wt') as h:
        for n,s in reads:h.write(f'@{n}\n{s}\n+\n{"I"*len(s)}\n')
def write_pe(p1,p2,kind):
    pairs=[]
    centers=[30000,80000] if kind=='chip' else []
    for i in range(1000):
        start=(centers[i%2]-400+(i*19)%800) if centers and i<750 else (i*883)%1985000; insert=180+(i%40)
        pairs.append((f'{kind}pe_{i}',seq[start:start+50],rc(seq[start+insert-50:start+insert])))
    for i in range(20):
        start=5100+i; insert=180
        pairs.append((f'{kind}pe_blacklist_{i}',seq[start:start+50],rc(seq[start+insert-50:start+insert])))
    for i in range(20):
        start=100+i; insert=180
        pairs.append((f'{kind}pe_mt_{i}',mt[start:start+50],rc(mt[start+insert-50:start+insert])))
    pairs.extend((f'{n}_duplicate_{i}',x,y) for i,(n,x,y) in enumerate(pairs[:15]))
    with gzip.open(p1,'wt') as a,gzip.open(p2,'wt') as b:
        for n,x,y in pairs:a.write(f'@{n}/1\n{x}\n+\n{"I"*50}\n');b.write(f'@{n}/2\n{y}\n+\n{"I"*50}\n')
write_se(raw/'se_chip.fastq.gz','chip');write_se(raw/'se_input.fastq.gz','input');write_pe(raw/'pe_chip_R1.fastq.gz',raw/'pe_chip_R2.fastq.gz','chip');write_pe(raw/'pe_input_R1.fastq.gz',raw/'pe_input_R2.fastq.gz','input')
cols=['sample_id','library_id','subject_id','donor_id','biosample_id','condition','assay','biological_replicate','technical_replicate','batch','sex','treatment','timepoint','genotype','library_layout','strandedness','bisulfite_directionality','bisulfite_protocol','conversion_control_contig','read1','read2','processed_input','processed_input_type','source_genome_build','source_annotation','source_processing_pipeline','source_file_accession','source_file_checksum','coordinate_normalization_status','gene_id_mapping_status','control_sample_id','histone_target','genome_build']
def blank():return {k:'' for k in cols}
def rawrow(s,l,assay,layout,r1,r2,control,target,rep):
    x=blank();x.update(sample_id=s,library_id=l,biosample_id='synthetic_biosample',condition='baseline',assay=assay,biological_replicate=str(rep),technical_replicate='1',batch='batch1',treatment='none',timepoint='0h',genotype='synthetic',library_layout=layout,read1=wslpath(r1),read2=wslpath(r2) if r2 else '',control_sample_id=control,histone_target=target,genome_build='GRCh38_primary');return x
for name,rows in [('raw_single.tsv',[rawrow('chip_se','chip_se_lib','H3K27ac_ChIP','single',raw/'se_chip.fastq.gz',None,'input_se','H3K27ac',1),rawrow('input_se','input_se_lib','ChIP_input','single',raw/'se_input.fastq.gz',None,'','CONTROL',1)]),('raw_paired.tsv',[rawrow('chip_pe','chip_pe_lib','H3K27ac_ChIP','paired',raw/'pe_chip_R1.fastq.gz',raw/'pe_chip_R2.fastq.gz','input_pe','H3K27ac',1),rawrow('input_pe','input_pe_lib','ChIP_input','paired',raw/'pe_input_R1.fastq.gz',raw/'pe_input_R2.fastq.gz','','CONTROL',1)])]:
    with open(root/name,'w',newline='') as h:w=csv.DictWriter(h,cols,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
samples=[]
for condition in ('control','treatment'):
 for i in range(1,4):
  sid=f'{condition}_{i}';cid=f'{sid}_input';samples.extend([(sid,'H3K27ac_ChIP','H3K27ac',cid),(cid,'ChIP_input','CONTROL','')])
regions=[('chr1',1000,1200),('chr1',3000,3200),('chr1',5000,5200),('chr1',7000,7200),('chr1',9000,9200)]
with open(proc/'counts.tsv','w',newline='') as h:
 w=csv.writer(h,delimiter='\t',lineterminator='\n');w.writerow(['region_id','chromosome','start','end']+[x[0] for x in samples])
 for j,(c,s,e) in enumerate(regions):
  rid='H3K27AC_'+hashlib.sha256(f'{c}:{s}-{e}'.encode()).hexdigest()[:16];vals=[]
  for sid,assay,tar,cid in samples:
   if assay=='ChIP_input':v=20+j
   elif sid.startswith('control'):v=[100,1000,600,2,0][j]+int(sid[-1])
   else:v=0 if j==4 else [100,8000,50,2,0][j]+int(sid[-1])
   if assay=='H3K27ac_ChIP' and j==4:v=0
   vals.append(v)
  w.writerow([rid,c,s,e]+vals)
rows=[]
for idx,(sid,assay,target,cid) in enumerate(samples):
 x=blank();condition=sid.split('_')[0];rep=sid.split('_')[1]
 x.update(sample_id=sid,library_id=sid+'_lib',biosample_id='synthetic_'+condition,condition=condition,assay=assay,biological_replicate=rep,technical_replicate='1',batch='batch1',treatment=condition,timepoint='0h',genotype='synthetic',library_layout='not_applicable',processed_input=wslpath(proc/'counts.tsv'),processed_input_type='h3k27ac_region_count_matrix_tsv',source_genome_build='synthetic_GRCh38_compatible',source_annotation='not_applicable',source_processing_pipeline='synthetic_counts_v1',source_file_accession='NOT_APPLICABLE',source_file_checksum='UNKNOWN',coordinate_normalization_status='verified_primary_assembly',gene_id_mapping_status='not_applicable',control_sample_id=cid,histone_target=target,genome_build='GRCh38_primary');rows.append(x)
with open(root/'processed.tsv','w',newline='') as h:w=csv.DictWriter(h,cols,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
with open(root/'one_vs_one.tsv','w',newline='') as h:
 w=csv.DictWriter(h,cols,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows([x for x in rows if x['sample_id'] in ('control_1','control_1_input','treatment_1','treatment_1_input')])
rank_rows=[]
for x in rows:
 y=dict(x);y['batch']=y['condition'];rank_rows.append(y)
with open(root/'rank_deficient.tsv','w',newline='') as h:w=csv.DictWriter(h,cols,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rank_rows)
