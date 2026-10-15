#!/usr/bin/env python3
import csv, json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]/'tests/output'
def status(path):
    with open(path) as h: return {r['component']:r['status'] for r in csv.DictReader(h,delimiter='\t')}
raw=ROOT/'atac_raw/atac'; proc=ROOT/'atac_processed/atac'
s=status(raw/'provenance/atac.analysis_status.tsv')
assert s['atac_processing']=='COMPLETE' and s['atac_peak_calling']=='COMPLETE' and s['atac_count_contract']=='COMPLETE' and s['atac_differential']=='NOT_RUN_DISABLED'
metrics={r['metric']:r['value'] for r in csv.DictReader(open(raw/'qc/alignment/atac_lib.filter_metrics.tsv'),delimiter='\t')}
assert int(metrics['all_first_mates'])==950 and int(metrics['mitochondrial'])==50 and int(metrics['duplicate'])>0 and float(metrics['NRF'])>0
qc=json.load(open(raw/'qc/metrics/atac_sample.atac_qc.json')); assert qc['blacklist_fragments']>0 and qc['usable_fragments']>0 and qc['peak_count']>=2 and qc['tss_enrichment'] is not None
regions=list(csv.DictReader(open(raw/'contract/atac.regions.tsv'),delimiter='\t')); assert len(regions)>=2 and len({r['region_id'] for r in regions})==len(regions)
matrix=list(csv.DictReader(open(raw/'counts/atac.raw_counts.tsv'),delimiter='\t')); assert all(r['atac_sample'].isdigit() for r in matrix) and sum(int(r['atac_sample']) for r in matrix)>0
assert (raw/'alignment/samples/atac_sample.bam').stat().st_size>0 and (raw/'alignment/samples/atac_sample.bam.bai').stat().st_size>0
assert (raw/'qc/ataqv/atac_sample.ataqv.json').stat().st_size>0 and (raw/'multiqc/atac_multiqc_report.html').stat().st_size>0 and (raw/'provenance/atac.provenance.json').stat().st_size>0
ps=status(proc/'provenance/atac.analysis_status.tsv'); assert ps['atac_count_contract']=='COMPLETE' and ps['atac_differential']=='COMPLETE'
de=list(csv.DictReader(open(proc/'differential/atac.differential.complete.tsv'),delimiter='\t')); assert len(de)==100
assert all(float(r['log2FoldChange'])>0 for r in de[:15]) and all(float(r['log2FoldChange'])<0 for r in de[15:30])
assert any(r['pvalue']=='NA' or r['padj']=='NA' for r in de)
print('ATAC output assertions passed')
