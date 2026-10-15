import csv,json,subprocess,sys,tempfile,unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]; FIX=ROOT/'tests'/'fixtures'/'generated'
class Contracts(unittest.TestCase):
    @classmethod
    def setUpClass(cls): subprocess.run([sys.executable,str(ROOT/'tests'/'create_synthetic_fixtures.py')],check=True,capture_output=True,text=True)
    def run_validator(self,name,ok=True):
        with tempfile.TemporaryDirectory() as d:
            p=subprocess.run([sys.executable,str(ROOT/'bin'/'validate_samplesheet.py'),'--input',str(FIX/name),'--output',str(Path(d)/'out.tsv'),'--report',str(Path(d)/'report.json')],capture_output=True,text=True,cwd=d)
            self.assertEqual(p.returncode==0,ok,p.stderr+p.stdout)
            if ok: return json.loads((Path(d)/'report.json').read_text())
    def test_raw_valid(self): self.assertEqual(self.run_validator('raw_samplesheet.tsv')['raw_libraries'],1)
    def test_processed_valid(self): self.assertEqual(self.run_validator('processed_samplesheet.tsv')['biological_samples'],6)
    def test_raw_processed_conflict(self): self.run_validator('raw_processed_conflict.tsv',False)
    def test_malformed_pair(self): self.run_validator('malformed_paired.tsv',False)
    def test_reference(self):
        with tempfile.TemporaryDirectory() as d:
            p=subprocess.run([sys.executable,str(ROOT/'bin'/'validate_reference.py'),'--fasta',str(FIX/'reference.fa'),'--gtf',str(FIX/'reference.gtf'),'--scope','synthetic_execution_only','--output',str(Path(d)/'r.json')],capture_output=True,text=True)
            self.assertEqual(p.returncode,0,p.stderr); self.assertEqual(json.loads((Path(d)/'r.json').read_text())['fasta']['contigs'],['chr22'])
if __name__=='__main__': unittest.main()
