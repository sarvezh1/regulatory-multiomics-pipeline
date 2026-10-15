#!/usr/bin/env python3
import csv, json
from pathlib import Path

OUT=Path(__file__).resolve().parents[1]/'tests/output'
def status(root):
    return {r['component']:r['status'] for r in csv.DictReader(open(root/'provenance/h3k27ac.analysis_status.tsv'),delimiter='\t')}

for name in ('chip_raw_single','chip_raw_paired'):
    root=OUT/name/'h3k27ac'; states=status(root)
    assert states['h3k27ac_processing']=='COMPLETE'
    assert states['h3k27ac_peak_calling']=='COMPLETE'
    assert states['h3k27ac_count_contract']=='COMPLETE'
    assert states['h3k27ac_differential']=='NOT_RUN_DISABLED'
    qc=json.load(open(next((root/'qc/metrics').glob('*.chip_qc.json'))))
    assert qc['mapped_reads_or_fragments']>0 and qc['duplicate']>0
    assert qc['mitochondrial']>0 and qc['blacklist_removed']>0 and qc['peak_count']>0
    regions=list(csv.DictReader(open(root/'contract/h3k27ac.regions.tsv'),delimiter='\t'))
    assert regions and len({r['region_id'] for r in regions})==len(regions)
    counts=list(csv.DictReader(open(root/'counts/h3k27ac.raw_counts.tsv'),delimiter='\t'))
    sample=[x for x in counts[0] if x not in ('region_id','chromosome','start','end')][0]
    assert all(r[sample].isdigit() for r in counts)
    assert (root/'multiqc/h3k27ac_multiqc_report.html').stat().st_size>0
    assert (root/'provenance/h3k27ac.provenance.json').stat().st_size>0

single=OUT/'chip_raw_single/h3k27ac'
assert 'predicted fragment length is' in (single/'peaks/samples/chip_se.macs3.log').read_text()

proc=OUT/'chip_processed/h3k27ac'; states=status(proc)
assert states['h3k27ac_count_contract']=='COMPLETE' and states['h3k27ac_differential']=='COMPLETE'
rows=list(csv.DictReader(open(proc/'differential/h3k27ac.differential.complete.tsv'),delimiter='\t'))
assert len(rows)==5
assert float(rows[1]['log2FoldChange'])>0 and float(rows[2]['log2FoldChange'])<0
assert rows[4]['pvalue']=='NA' and rows[4]['padj']=='NA'
print('H3K27ac output assertions passed')
