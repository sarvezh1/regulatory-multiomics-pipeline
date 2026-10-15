import csv, subprocess, sys, tempfile, unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]; FIX=ROOT/'tests'/'fixtures'/'wgbs'
class WgbsContracts(unittest.TestCase):
    @classmethod
    def setUpClass(cls): subprocess.run([sys.executable,str(ROOT/'tests'/'create_wgbs_fixtures.py')],check=True)
    def validate(self,name,ok):
        with tempfile.TemporaryDirectory() as tmp:
            run=subprocess.run([sys.executable,str(ROOT/'bin'/'validate_wgbs_samplesheet.py'),'--input',str(FIX/name),'--output',str(Path(tmp)/'out.tsv'),'--report',str(Path(tmp)/'report.json')],capture_output=True,text=True)
            self.assertEqual(run.returncode==0,ok,run.stderr+run.stdout)
    def test_raw_valid(self): self.validate('raw_samplesheet.tsv',True)
    def test_processed_valid(self): self.validate('processed_samplesheet.tsv',True)
    def test_missing_r2(self): self.validate('missing_r2.tsv',False)
    def test_unknown_directionality(self): self.validate('unknown_directionality.tsv',False)
    def test_unsupported_protocol(self): self.validate('unsupported_protocol.tsv',False)
    def test_incompatible_genome(self): self.validate('incompatible_genome.tsv',False)
    def test_mixed_processed_types(self): self.validate('mixed_processed_types.tsv',False)
    def adapter(self, values, ok):
        fields=['chromosome','start','end','coordinate_system','methylated_count','unmethylated_count','coverage','methylation_fraction']
        with tempfile.TemporaryDirectory() as tmp:
            source=Path(tmp)/'input.tsv'; output=Path(tmp)/'output.tsv'
            with source.open('w',newline='',encoding='utf-8') as handle:
                writer=csv.DictWriter(handle,fieldnames=fields,delimiter='\t',lineterminator='\n'); writer.writeheader(); writer.writerow(dict(zip(fields,values)))
            run=subprocess.run([sys.executable,str(ROOT/'bin'/'adapt_processed_wgbs.py'),'--input',str(source),'--fasta',str(FIX/'processed_reference.fa'),'--sample-id','s','--library-id','l','--type','cpg_count_tsv','--output',str(output)],capture_output=True,text=True)
            self.assertEqual(run.returncode==0,ok,run.stderr+run.stdout)
    def test_negative_count_rejected(self): self.adapter(['chr22',100,102,'0-based-half-open',-1,11,10,'NA'],False)
    def test_methylated_exceeds_coverage_rejected(self): self.adapter(['chr22',100,102,'0-based-half-open',11,'',10,'NA'],False)
    def test_fraction_outside_range_rejected(self): self.adapter(['chr22',100,102,'0-based-half-open',5,5,10,1.2],False)
    def test_invalid_coordinates_rejected(self): self.adapter(['chr22',100,100,'0-based-half-open',5,5,10,0.5],False)
    def test_zero_coverage_fraction_rejected(self): self.adapter(['chr22',100,102,'0-based-half-open',0,0,0,0],False)
    def test_coordinate_and_strand_merge(self):
        with tempfile.TemporaryDirectory() as tmp:
            cov=Path(tmp)/'x.cov'; cov.write_text('chr22\t101\t101\t50\t2\t2\nchr22\t102\t102\t50\t3\t3\n')
            out=Path(tmp)/'out.tsv'; rejected=Path(tmp)/'rejected.tsv'; subprocess.run([sys.executable,str(ROOT/'bin'/'wgbs_cpg_contract.py'),'--coverage',str(cov),'--fasta',str(FIX/'processed_reference.fa'),'--sample-id','s','--library-id','l','--output',str(out),'--rejected-output',str(rejected)],check=True)
            with out.open() as handle:
                rows=list(csv.DictReader(handle,delimiter='\t'))
            self.assertEqual(len(rows),1); self.assertEqual((rows[0]['start'],rows[0]['end']),('100','102')); self.assertEqual(rows[0]['coverage'],'10'); self.assertEqual(rows[0]['methylated_count'],'5')
            with rejected.open() as handle: excluded=list(csv.DictReader(handle,delimiter='\t'))
            self.assertEqual(len(excluded),0)
    def test_nonreference_cpg_is_explicitly_excluded(self):
        with tempfile.TemporaryDirectory() as tmp:
            cov=Path(tmp)/'x.cov'; cov.write_text('chr22\t105\t105\t0\t0\t1\n')
            out=Path(tmp)/'out.tsv'; rejected=Path(tmp)/'rejected.tsv'
            subprocess.run([sys.executable,str(ROOT/'bin'/'wgbs_cpg_contract.py'),'--coverage',str(cov),'--fasta',str(FIX/'processed_reference.fa'),'--sample-id','s','--library-id','l','--output',str(out),'--rejected-output',str(rejected)],check=True)
            with rejected.open() as handle: excluded=list(csv.DictReader(handle,delimiter='\t'))
            self.assertEqual(len(excluded),1); self.assertEqual(excluded[0]['reason'].split(':')[0],'non-CpG source coordinate')

if __name__=='__main__': unittest.main()
