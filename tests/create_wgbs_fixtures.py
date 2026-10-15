#!/usr/bin/env python3
import csv, gzip, hashlib, random
from pathlib import Path

ROOT=Path(__file__).resolve().parent/'fixtures'/'wgbs'
HEADERS=['sample_id','library_id','subject_id','donor_id','biosample_id','condition','assay','biological_replicate','technical_replicate','batch','sex','treatment','timepoint','genotype','library_layout','strandedness','bisulfite_directionality','bisulfite_protocol','conversion_control_contig','read1','read2','processed_input','processed_input_type','source_genome_build','source_annotation','source_processing_pipeline','source_file_accession','source_file_checksum','coordinate_normalization_status','gene_id_mapping_status','control_sample_id','histone_target','genome_build']

def revcomp(seq): return seq.translate(str.maketrans('ACGT','TGCA'))[::-1]
def converted(seq,methylated_positions,offset):
    chars=list(seq)
    for i,base in enumerate(chars):
        if base=='C' and offset+i not in methylated_positions: chars[i]='T'
    return ''.join(chars)
def converted_by_index(seq,methylated_indices):
    chars=list(seq)
    for i,base in enumerate(chars):
        if base=='C' and i not in methylated_indices: chars[i]='T'
    return ''.join(chars)
def fq(path,records):
    with gzip.GzipFile(path,'wb',mtime=0) as gz:
        for name,seq in records: gz.write(f'@{name}\n{seq}\n+\n{"I"*len(seq)}\n'.encode())
def base_row(**kwargs):
    row={h:'' for h in HEADERS}; row.update({'subject_id':'synthetic_subject','donor_id':'synthetic_donor','biosample_id':'synthetic_biosample','assay':'WGBS','batch':'batch1','sex':'unknown','treatment':'none','timepoint':'0 h','genotype':'parent','strandedness':'','bisulfite_protocol':'conventional_wgbs','source_annotation':'NOT_APPLICABLE','gene_id_mapping_status':'not_applicable','control_sample_id':'','histone_target':'','genome_build':'GRCh38_primary'}); row.update(kwargs); return row
def write_sheet(name,rows):
    with (ROOT/name).open('w',newline='',encoding='utf-8') as out:
        w=csv.DictWriter(out,fieldnames=HEADERS,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(rows)

def main():
    ROOT.mkdir(parents=True,exist_ok=True); (ROOT/'raw').mkdir(exist_ok=True); (ROOT/'processed').mkdir(exist_ok=True)
    rng=random.Random(31003); genome=[rng.choice('ACGT') for _ in range(6000)]
    # Force known CpGs and keep unique flanks.
    cpgs=[900,930,960,1200,1230,1260,2200,2230,2260]
    for p in cpgs: genome[p:p+2]=list('CG')
    reference=''.join(genome)
    lambda_seq=('ATGCGTACGTTAGCTACGATCGATGCTAGCATCGATCGTACGATCGAT'*20)[:900]
    with (ROOT/'reference.fa').open('w',encoding='utf-8') as out: out.write('>chr22\n'+reference+'\n>chrM\n'+lambda_seq+'\n')
    with (ROOT/'reference.gtf').open('w',encoding='utf-8') as out: out.write('chr22\tfixture\tgene\t801\t1400\t.\t+\t.\tgene_id "GENE1.1"; gene_name "GENE1";\nchr22\tfixture\ttranscript\t801\t1400\t.\t+\t.\tgene_id "GENE1.1"; transcript_id "TX1.1"; gene_name "GENE1";\nchr22\tfixture\texon\t801\t1400\t.\t+\t.\tgene_id "GENE1.1"; transcript_id "TX1.1"; gene_name "GENE1";\n')
    raw_rows=[]
    for lib,tech,starts in [('wgbs_lib1','1',[850,880,910,1150,1180,1210,2150,2180]),('wgbs_lib2','2',[860,890,920,1160,1190,1220,2160,2190])]:
        r1=[]; r2=[]; idx=0
        for start in starts:
            for rep in range(4):
                insert=reference[start:start+180]; local_cpg={p-start for p in cpgs if start<=p<start+180}
                if rep < 2:
                    methylated={start+p for p in local_cpg if ((p//30)+rep)%2==0}
                    molecule=converted(insert,methylated,start)
                else:
                    bottom=revcomp(insert)
                    methylated_bottom={len(insert)-1-(p+1-start) for p in cpgs if start<=p<start+179 and ((p//30)+rep)%2==0}
                    molecule=converted_by_index(bottom,methylated_bottom)
                name=f'{lib}_{idx:04d}'; r1.append((name+'/1',molecule[:80])); r2.append((name+'/2',revcomp(molecule[-80:]))); idx+=1
        # Two PCR duplicate pairs share genomic coordinates and sequence with the first pair.
        r1.extend([(f'{lib}_dup{i}/1',r1[0][1]) for i in range(2)]); r2.extend([(f'{lib}_dup{i}/2',r2[0][1]) for i in range(2)])
        fq(ROOT/'raw'/f'{lib}_R1.fastq.gz',r1); fq(ROOT/'raw'/f'{lib}_R2.fastq.gz',r2)
        raw_rows.append(base_row(sample_id='wgbs_sample',library_id=lib,condition='control',biological_replicate='1',technical_replicate=tech,library_layout='paired',bisulfite_directionality='directional',conversion_control_contig='',read1=f'raw/{lib}_R1.fastq.gz',read2=f'raw/{lib}_R2.fastq.gz'))
    write_sheet('raw_samplesheet.tsv',raw_rows)
    fq(ROOT/'raw'/'zero_align_R1.fastq.gz',[('zero/1','A'*80)]); fq(ROOT/'raw'/'zero_align_R2.fastq.gz',[('zero/2','A'*80)])
    zero_raw=base_row(sample_id='zero_align',library_id='zero_align_lib',condition='control',biological_replicate='1',technical_replicate='1',library_layout='paired',bisulfite_directionality='directional',conversion_control_contig='',read1='raw/zero_align_R1.fastq.gz',read2='raw/zero_align_R2.fastq.gz')
    write_sheet('zero_align.tsv',[zero_raw])

    samples=[]; positions=list(range(100,1100,10))
    for group in ('control','treatment'):
        for rep in range(1,4):
            sid=f'{group}_{rep}'; path=ROOT/'processed'/f'{sid}.tsv'
            with path.open('w',newline='',encoding='utf-8') as out:
                fields=['chromosome','start','end','coordinate_system','methylated_count','unmethylated_count','coverage','methylation_fraction']; w=csv.DictWriter(out,fieldnames=fields,delimiter='\t',lineterminator='\n'); w.writeheader()
                for i,pos in enumerate(positions):
                    # Reference-aware adapter uses chrSynthetic; insert CpGs at these positions in a dedicated processed reference below.
                    cov=10
                    if 10<=i<35: meth=2 if group=='control' else 9
                    elif 50<=i<75: meth=8 if group=='control' else 1
                    else: meth=4+(rep%2)
                    if i in (80,81): cov=2; meth=min(meth,cov)
                    if i==90 and rep==3: continue
                    fraction=f'{meth/cov:.8f}'
                    w.writerow({'chromosome':'chr22','start':pos,'end':pos+2,'coordinate_system':'0-based-half-open','methylated_count':meth,'unmethylated_count':cov-meth,'coverage':cov,'methylation_fraction':fraction})
            samples.append(base_row(sample_id=sid,library_id=sid+'_lib',condition=group,biological_replicate=str(rep),technical_replicate='1',library_layout='not_applicable',bisulfite_directionality='UNKNOWN',processed_input=f'processed/{sid}.tsv',processed_input_type='cpg_count_tsv',source_genome_build='synthetic_GRCh38_compatible',source_processing_pipeline='deterministic_fixture',source_file_accession='NOT_APPLICABLE',source_file_checksum='UNKNOWN',coordinate_normalization_status='verified_primary_assembly'))
    # Give every processed locus a valid reference CpG.
    genome2=list(reference)
    for pos in positions: genome2[pos:pos+2]=list('CG')
    with (ROOT/'processed_reference.fa').open('w',encoding='utf-8') as out: out.write('>chr22\n'+''.join(genome2)+'\n')
    write_sheet('processed_samplesheet.tsv',samples)
    write_sheet('one_vs_one.tsv',[samples[0],samples[3]])
    rank=[]
    for row in samples:
        item=dict(row); item['batch']=item['condition']; rank.append(item)
    write_sheet('rank_deficient.tsv',rank)
    for mode,cov in [('low_coverage',1),('zero_coverage',0)]:
        mode_rows=[]
        for group in ('control','treatment'):
            for rep in (1,2):
                sid=f'{mode}_{group}_{rep}'; path=ROOT/'processed'/f'{sid}.tsv'
                with path.open('w',newline='',encoding='utf-8') as out:
                    fields=['chromosome','start','end','coordinate_system','methylated_count','unmethylated_count','coverage','methylation_fraction']; w=csv.DictWriter(out,fieldnames=fields,delimiter='\t',lineterminator='\n'); w.writeheader()
                    for pos in positions[:10]: w.writerow({'chromosome':'chr22','start':pos,'end':pos+2,'coordinate_system':'0-based-half-open','methylated_count':0,'unmethylated_count':cov,'coverage':cov,'methylation_fraction':'0.00000000' if cov else 'NA'})
                mode_rows.append(base_row(sample_id=sid,library_id=sid+'_lib',condition=group,biological_replicate=str(rep),technical_replicate='1',library_layout='not_applicable',bisulfite_directionality='UNKNOWN',processed_input=f'processed/{sid}.tsv',processed_input_type='cpg_count_tsv',source_genome_build='synthetic_GRCh38_compatible',source_processing_pipeline='deterministic_fixture',source_file_accession='NOT_APPLICABLE',source_file_checksum='UNKNOWN',coordinate_normalization_status='verified_primary_assembly'))
        write_sheet(mode+'.tsv',mode_rows)
    # Percentage-only descriptive route.
    percent=ROOT/'processed'/'percentage.tsv'
    with percent.open('w',newline='',encoding='utf-8') as out:
        w=csv.DictWriter(out,fieldnames=['chromosome','start','end','coordinate_system','methylation_fraction'],delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerow({'chromosome':'chr22','start':100,'end':102,'coordinate_system':'0-based-half-open','methylation_fraction':'0.5'})
    p_rows=[]
    for group,rep in [('control',1),('control',2),('treatment',1),('treatment',2)]: p_rows.append(base_row(sample_id=f'pct_{group}_{rep}',library_id=f'pct_{group}_{rep}_lib',condition=group,biological_replicate=str(rep),technical_replicate='1',library_layout='not_applicable',bisulfite_directionality='UNKNOWN',processed_input='processed/percentage.tsv',processed_input_type='cpg_percentage_tsv',source_genome_build='synthetic_GRCh38_compatible',source_processing_pipeline='external_percentage_fixture',source_file_accession='NOT_APPLICABLE',source_file_checksum='UNKNOWN',coordinate_normalization_status='verified_primary_assembly'))
    write_sheet('percentage_only.tsv',p_rows)
    incompatible=dict(samples[0]); incompatible['source_genome_build']='GRCm39'; incompatible['coordinate_normalization_status']='UNRESOLVED'; write_sheet('incompatible_genome.tsv',[incompatible])
    mixed=[dict(samples[0]),dict(p_rows[0])]; write_sheet('mixed_processed_types.tsv',mixed)
    bad=dict(raw_rows[0]); bad['read2']=''; write_sheet('missing_r2.tsv',[bad])
    bad=dict(raw_rows[0]); bad['bisulfite_directionality']='UNKNOWN'; write_sheet('unknown_directionality.tsv',[bad])
    bad=dict(raw_rows[0]); bad['bisulfite_protocol']='RRBS'; write_sheet('unsupported_protocol.tsv',[bad])

if __name__=='__main__': main()
