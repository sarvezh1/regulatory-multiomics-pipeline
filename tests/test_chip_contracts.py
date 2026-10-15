import csv, json, subprocess, sys, tempfile, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FIX = ROOT / 'tests/fixtures/chip'

def run(*args):
    return subprocess.run([sys.executable, *map(str, args)], cwd=ROOT, text=True, capture_output=True)

class ChipContracts(unittest.TestCase):
    def validate(self, source, ok=True):
        with tempfile.TemporaryDirectory() as d:
            result = run('bin/validate_chip_samplesheet.py', '--input', source,
                         '--output', Path(d)/'valid.tsv', '--report', Path(d)/'report.json')
            self.assertEqual(result.returncode == 0, ok, result.stderr + result.stdout)

    def mutated(self, source, mutate):
        with open(source) as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        for row in rows:
            for field in ('read1','read2','processed_input'):
                value=row.get(field,'')
                if value.startswith('/mnt/d/'):
                    row[field]='D:/'+value[len('/mnt/d/'):]
        mutate(rows)
        handle = tempfile.NamedTemporaryFile('w', suffix='.tsv', newline='', delete=False)
        writer = csv.DictWriter(handle, rows[0].keys(), delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(rows); handle.close()
        self.addCleanup(Path(handle.name).unlink)
        return Path(handle.name)

    def test_valid_routes(self):
        for name in ('raw_single.tsv', 'raw_paired.tsv', 'processed.tsv'):
            self.validate(self.mutated(FIX/name, lambda rows: None))

    def test_control_and_target_guards(self):
        base = FIX/'raw_single.tsv'
        mutations = [
            lambda r: r[0].update(control_sample_id=''),
            lambda r: r[0].update(control_sample_id='absent'),
            lambda r: r[0].update(histone_target='H3K4me3'),
            lambda r: r[1].update(condition='other'),
        ]
        for mutation in mutations:
            self.validate(self.mutated(base, mutation), False)

    def test_layout_and_input_guards(self):
        paired = FIX/'raw_paired.tsv'; single = FIX/'raw_single.tsv'
        self.validate(self.mutated(paired, lambda r: r[0].update(read2='')), False)
        self.validate(self.mutated(single, lambda r: r[0].update(read2=r[0]['read1'])), False)
        self.validate(self.mutated(single, lambda r: r[0].update(processed_input=r[0]['read1'], processed_input_type='h3k27ac_region_count_matrix_tsv')), False)
        self.validate(self.mutated(single, lambda r: r[0].update(genome_build='hg19')), False)

    def test_zero_peak_status(self):
        with tempfile.TemporaryDirectory() as d:
            empty=Path(d)/'empty.narrowPeak'; empty.write_text('')
            output=Path(d)/'status.json'
            result=run('bin/build_chip_status.py','--route','raw','--peaks',empty,'--output',output)
            self.assertEqual(result.returncode,0,result.stderr)
            states={x['component']:x['status'] for x in json.loads(output.read_text())['components']}
            self.assertEqual(states['h3k27ac_peak_calling'],'COMPLETE_EMPTY_RESULT')

    def test_consensus_strict_overlap_and_stable_ids(self):
        with tempfile.TemporaryDirectory() as d:
            peaks=Path(d)/'sample_peaks.narrowPeak'
            peaks.write_text('chr1\t10\t20\ta\t1\t.\t1\t1\t1\t5\nchr1\t19\t30\tb\t1\t.\t1\t1\t1\t5\nchr1\t30\t40\tc\t1\t.\t1\t1\t1\t5\n')
            regions=Path(d)/'regions.tsv'; provenance=Path(d)/'provenance.tsv'
            result=run('bin/build_chip_consensus.py','--peaks',peaks,'--regions',regions,'--provenance',provenance)
            self.assertEqual(result.returncode,0,result.stderr)
            with open(regions) as handle:
                rows=list(csv.DictReader(handle,delimiter='\t'))
            self.assertEqual([(x['start'],x['end']) for x in rows],[('10','30'),('30','40')])
            self.assertEqual(len({x['region_id'] for x in rows}),2)

if __name__ == '__main__': unittest.main()
