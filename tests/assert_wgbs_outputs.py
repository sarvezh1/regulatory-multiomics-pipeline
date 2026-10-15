#!/usr/bin/env python3
import csv, json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]/'tests'/'output'

def rows(path):
    with path.open(encoding='utf-8') as handle:
        return list(csv.DictReader(handle,delimiter='\t'))

def status(root):
    return {row['component']:row['status'] for row in rows(root/'provenance'/'wgbs.analysis_status.tsv')}

raw=ROOT/'wgbs_raw'/'wgbs'
raw_status=status(raw)
assert raw_status=={
    'wgbs_processing':'COMPLETE',
    'wgbs_methylation_extraction':'COMPLETE',
    'wgbs_count_contract':'COMPLETE',
    'wgbs_differential':'NOT_RUN_DISABLED',
}
qc=[]
for path in sorted((raw/'qc'/'metrics').glob('*.wgbs_qc.tsv')):
    qc.extend(rows(path))
assert len(qc)==2
assert all(int(row['input_pairs'])==34 and int(row['retained_pairs'])==34 for row in qc)
assert all(int(row['unique_aligned_pairs'])>0 and int(row['duplicates_removed'])>0 for row in qc)
assert all(row['mbias_status']=='COMPLETE' and row['conversion_control_status']=='NOT_AVAILABLE_INPUT' for row in qc)
sample_rows=rows(raw/'counts'/'samples'/'wgbs_sample.cpg.tsv')
assert sample_rows and all(int(row['methylated_count'])+int(row['unmethylated_count'])==int(row['coverage']) for row in sample_rows)
assert all(row['coordinate_system']=='0-based-half-open' and int(row['end'])-int(row['start'])==2 for row in sample_rows)
assert any(row['source_strands']=='bottom,top' for row in sample_rows)
assert all(row['contributing_libraries']=='wgbs_lib1,wgbs_lib2' for row in sample_rows)
assert (raw/'multiqc'/'wgbs_multiqc_report.html').stat().st_size>0
assert (raw/'provenance'/'wgbs.provenance.json').stat().st_size>0

processed=ROOT/'wgbs_processed'/'wgbs'
processed_status=status(processed)
assert processed_status['wgbs_count_contract']=='COMPLETE'
assert processed_status['wgbs_differential']=='COMPLETE'
dmr=rows(processed/'differential'/'wgbs.dmrseq.complete.tsv')
assert dmr and len({row['dmr_id'] for row in dmr})==len(dmr)
assert any(float(row['delta_methylation'])>0.5 for row in dmr)
assert any(float(row['delta_methylation'])<-0.5 for row in dmr)
long_rows=rows(processed/'counts'/'wgbs.cpg.long.tsv')
assert any(row['measurement_state']=='MISSING' for row in long_rows)

expected={
    'wgbs_guard_1v1':'NOT_RUN_INSUFFICIENT_REPLICATION',
    'wgbs_guard_no_contrast':'NOT_RUN_NO_CONTRAST',
    'wgbs_guard_rank':'NOT_RUN_INVALID_DESIGN',
    'wgbs_guard_percentage':'NOT_AVAILABLE_INPUT',
    'wgbs_guard_low_coverage':'COMPLETE_EMPTY_RESULT',
    'wgbs_guard_zero_coverage':'COMPLETE_EMPTY_RESULT',
}
for directory,state in expected.items():
    assert status(ROOT/directory/'wgbs')['wgbs_differential']==state

zero_raw=status(ROOT/'wgbs_zero_align'/'wgbs')
assert zero_raw['wgbs_processing']=='COMPLETE'
assert zero_raw['wgbs_methylation_extraction']=='COMPLETE_EMPTY_RESULT'

genuine=ROOT/'wgbs_genuine'/'wgbs'
genuine_status=status(genuine)
assert genuine_status=={
    'wgbs_processing':'COMPLETE',
    'wgbs_methylation_extraction':'COMPLETE',
    'wgbs_count_contract':'COMPLETE',
    'wgbs_differential':'NOT_RUN_DISABLED',
}
genuine_qc=rows(genuine/'qc'/'metrics'/'SRR4235744.wgbs_qc.tsv')
assert len(genuine_qc)==1
assert int(genuine_qc[0]['input_pairs'])==250000
assert int(genuine_qc[0]['retained_pairs'])==189997
assert int(genuine_qc[0]['unique_aligned_pairs'])==3303
genuine_provenance=json.loads((genuine/'provenance'/'wgbs.provenance.json').read_text())
assert genuine_provenance['reference']['reference_scope']=='compact_chr22_execution'
assert genuine_provenance['reference']['annotation_release']=='GENCODE 50'

print('WGBS output assertions passed')
