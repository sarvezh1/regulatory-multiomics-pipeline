#!/usr/bin/env python3
import argparse,csv,json
from pathlib import Path

def rows(path):
    with open(path,encoding='utf-8',newline='') as h: return list(csv.DictReader(h,delimiter='\t'))
def status(root,name): return json.loads((root/name/'rna/differential/rna.differential.status.json').read_text())['status']
def main():
    p=argparse.ArgumentParser(); p.add_argument('--root',default='tests/output'); p.add_argument('--include-genuine',action='store_true'); a=p.parse_args(); root=Path(a.root)
    raw=root/'raw'/'rna'; assert (raw/'alignment/raw_lib.Aligned.sortedByCoord.out.bam').stat().st_size>0; assert (raw/'alignment/raw_lib.Aligned.sortedByCoord.out.bam.bai').stat().st_size>0
    assert int((raw/'alignment/raw_lib.flagstat.txt').read_text().split()[0])>0
    matrix=rows(raw/'counts/rna.raw_counts.tsv'); assert {r['gene_name']:int(r['raw_sample']) for r in matrix}=={'UPTEST':12,'DOWNTEST':12}
    assert status(root,'raw')=='NOT_RUN_DISABLED'; assert (raw/'multiqc/rna_multiqc_report.html').stat().st_size>0; assert (raw/'provenance/rna.provenance.json').stat().st_size>0
    de=rows(root/'processed/rna/differential/rna.differential.complete.tsv'); by={r['gene_name']:r for r in de}; assert len(de)==100; assert float(by['GENE2']['log2FoldChange'])>0; assert float(by['GENE11']['log2FoldChange'])<0; assert by['GENE100']['pvalue']=='NA'
    assert status(root,'processed')=='COMPLETE'; assert status(root,'guard_1v1')=='NOT_RUN_INSUFFICIENT_REPLICATION'; assert status(root,'guard_no_contrast')=='NOT_RUN_NO_CONTRAST'; assert status(root,'guard_rank')=='NOT_RUN_INVALID_DESIGN'
    tech=rows(root/'technical_replicates/rna/counts/rna.technical_replicate_manifest.tsv'); meta=rows(root/'technical_replicates/rna/counts/rna.sample_metadata.tsv'); assert len(tech)==7 and len(meta)==6; assert status(root,'technical_replicates') in {'COMPLETE','COMPLETE_EMPTY_RESULT'}
    if a.include_genuine:
        genuine=root/'genuine'/'rna'; assert int((genuine/'alignment/SRR12413947_subset.flagstat.txt').read_text().split()[0])>0
        assert sum(int(r['GSM4715082'])>0 for r in rows(genuine/'counts/rna.raw_counts.tsv'))>0; assert status(root,'genuine')=='NOT_RUN_DISABLED'
    print('Output assertions passed')
if __name__=='__main__': main()
