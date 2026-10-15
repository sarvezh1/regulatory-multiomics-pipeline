import csv, json, subprocess, sys, tempfile, unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]; FIX=ROOT/'tests/fixtures/atac'
def run(*args): return subprocess.run([sys.executable,*map(str,args)],cwd=ROOT,text=True,capture_output=True)

class AtacContracts(unittest.TestCase):
    def validate(self,name,ok):
        with tempfile.TemporaryDirectory() as d:
            r=run('bin/validate_atac_samplesheet.py','--input',FIX/name,'--output',Path(d)/'x.tsv','--report',Path(d)/'x.json')
            self.assertEqual(r.returncode==0,ok,r.stderr+r.stdout)
    def test_raw_paired_valid(self): self.validate('raw_samplesheet.tsv',True)
    def test_processed_matrix_valid(self): self.validate('processed_samplesheet.tsv',True)
    def test_missing_r2(self): self.validate('missing_r2.tsv',False)
    def test_raw_processed_conflict(self): self.validate('raw_processed_conflict.tsv',False)
    def test_invalid_assay(self): self.validate('invalid_assay.tsv',False)
    def test_incompatible_build(self): self.validate('incompatible_genome.tsv',False)
    def test_malformed_bed(self):
        with tempfile.TemporaryDirectory() as d:
            frag=Path(d)/'f.bed'; frag.write_text('chr1\t10\t20\tq1\n')
            r=run('bin/filter_atac_fragments.py','--fragments',frag,'--blacklist',FIX/'malformed_blacklist.bed','--output',Path(d)/'o.bed','--names',Path(d)/'n.txt','--metrics',Path(d)/'m.json')
            self.assertNotEqual(r.returncode,0)
    def test_processed_technical_libraries_not_samples(self):
        with tempfile.TemporaryDirectory() as d:
            valid=Path(d)/'v.tsv'; report=Path(d)/'r.json'
            r=run('bin/validate_atac_samplesheet.py','--input',FIX/'processed_technical_duplicate.tsv','--output',valid,'--report',report); self.assertEqual(r.returncode,0,r.stderr)
            r=run('bin/adapt_atac_matrix.py','--matrix',FIX/'processed_atac_counts.tsv','--metadata',valid,'--output',Path(d)/'c.tsv','--regions',Path(d)/'g.tsv','--provenance',Path(d)/'p.json')
            self.assertNotEqual(r.returncode,0); self.assertIn('one row per biological sample',r.stderr+r.stdout)
    def test_zero_peak_status(self):
        with tempfile.TemporaryDirectory() as d:
            empty=Path(d)/'empty.narrowPeak'; empty.write_text('')
            output=Path(d)/'status.json'
            r=run('bin/build_atac_status.py','--route','raw','--peaks',empty,'--output',output)
            self.assertEqual(r.returncode,0,r.stderr)
            states={x['component']:x['status'] for x in json.loads(output.read_text())['components']}
            self.assertEqual(states['atac_peak_calling'],'COMPLETE_EMPTY_RESULT')
    def test_consensus_merges_overlap_not_book_ends(self):
        with tempfile.TemporaryDirectory() as d:
            peaks=Path(d)/'sample.narrowPeak'
            peaks.write_text('chr1\t10\t20\ta\nchr1\t19\t30\tb\nchr1\t30\t40\tc\n')
            regions=Path(d)/'regions.tsv'; provenance=Path(d)/'provenance.tsv'
            r=run('bin/build_atac_consensus.py','--peaks',peaks,'--regions',regions,'--provenance',provenance)
            self.assertEqual(r.returncode,0,r.stderr)
            with regions.open() as handle:
                rows=list(csv.DictReader(handle,delimiter='\t'))
            self.assertEqual([(x['start'],x['end']) for x in rows],[('10','30'),('30','40')])

if __name__=='__main__': unittest.main()
