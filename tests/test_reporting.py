import csv
import hashlib
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BUILD = ROOT / "bin" / "build_run_report.py"
VALIDATE = ROOT / "bin" / "validate_run_report.py"
SPEC = importlib.util.spec_from_file_location("build_run_report", BUILD)
REPORT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(REPORT)


def write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


class ReportingTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.source = self.root / "source"
        self.results = self.root / "results"
        self.report = self.root / "report"
        write(self.source / "main.nf", "nextflow.enable.dsl = 2\n")
        write(self.source / "nextflow.config", "params.project_version = 'test'\n")
        write(self.source / "assets" / "static.txt", "fixture\n")
        write(self.source / "sample.tsv", "sample_id\tlibrary_id\tassay\tread1\tsource_genome_build\nS1\tL1\tRNA\tNA\tGRCh38\n")
        write(self.source / "reference.fa", ">chr22\nACGT\n")
        write(self.source / "parameters.json", '{"run_de": false, "strandedness": "unstranded"}\n')
        digest = "sha256:" + "1" * 64
        write(self.source / "tools.tsv", f"process\ttool\tcontainer\trepository_digest\nTEST\tTool 1.0\texample/tool:1.0@{digest}\t{digest}\n")
        write(self.results / "rna" / "counts" / "rna.raw_counts.tsv", "gene_id\tS1\nENSG000001\t5\n")
        write(self.results / "rna" / "provenance" / "rna.analysis_status.tsv", "component\tstatus\nrna_count_contract\tCOMPLETE\nrna_differential\tNOT_RUN_DISABLED\n")

    def tearDown(self):
        self.temp.cleanup()

    def command(self, output=None):
        return [sys.executable, str(BUILD), "--results", str(self.results), "--source-root", str(self.source), "--output", str(output or self.report), "--requested-assay", "RNA", "--nextflow-version", "26.04.6", "--input-record", str(self.source / "sample.tsv"), "--reference", f"FASTA={self.source / 'reference.fa'}", "--parameter-file", str(self.source / "parameters.json"), "--tool-containers", str(self.source / "tools.tsv"), "--timestamp-utc", "2026-01-01T00:00:00+00:00"]

    def build(self):
        subprocess.run(self.command(), check=True, capture_output=True, text=True)

    def test_valid_report_and_optional_disabled_status(self):
        self.build()
        completed = subprocess.run([sys.executable, str(VALIDATE), str(self.report)], check=True, capture_output=True, text=True)
        self.assertIn("validation passed", completed.stdout)
        statuses = (self.report / "status_manifest.tsv").read_text(encoding="utf-8")
        self.assertIn("NOT_RUN_DISABLED", statuses)

    def test_missing_completed_output_is_rejected(self):
        (self.results / "rna" / "counts" / "rna.raw_counts.tsv").unlink()
        completed = subprocess.run(self.command(), capture_output=True, text=True)
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("required output is missing", completed.stderr)

    def test_malformed_tsv_and_json_are_rejected(self):
        write(self.source / "sample.tsv", "sample_id\tlibrary_id\nS1\n")
        completed = subprocess.run(self.command(), capture_output=True, text=True)
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("malformed TSV", completed.stderr)
        write(self.source / "sample.tsv", "sample_id\tlibrary_id\tassay\tread1\nS1\tL1\tRNA\tNA\n")
        write(self.source / "parameters.json", "{not-json}\n")
        completed = subprocess.run(self.command(), capture_output=True, text=True)
        self.assertNotEqual(completed.returncode, 0)

    def test_checksum_mismatch_is_rejected(self):
        self.build()
        with (self.report / "status_manifest.tsv").open("a", encoding="utf-8") as handle:
            handle.write("status_manifest\t1.0.0\tx\tCOMPLETE\tx\tx\n")
        completed = subprocess.run([sys.executable, str(VALIDATE), str(self.report)], capture_output=True, text=True)
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("checksum mismatch", completed.stderr)

    def test_private_path_is_rejected(self):
        private_value = "C:" + "\\\\Users\\\\person\\\\reads.fastq.gz"
        write(self.source / "parameters.json", json.dumps({"private_input": private_value}) + "\n")
        completed = subprocess.run(self.command(), capture_output=True, text=True)
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("private path", completed.stderr)

    def test_source_fingerprint_changes_with_source(self):
        _, before = REPORT.source_manifest(self.source)
        write(self.source / "bin" / "new_step.py", "print('changed')\n")
        _, after = REPORT.source_manifest(self.source)
        self.assertNotEqual(before, after)

    def test_duplicate_output_identifier_is_rejected(self):
        self.build()
        output_path = self.report / "output_manifest.tsv"
        lines = output_path.read_text(encoding="utf-8").splitlines()
        output_path.write_text("\n".join(lines + [lines[1]]) + "\n", encoding="utf-8")
        run_path = self.report / "run_manifest.json"
        run = json.loads(run_path.read_text(encoding="utf-8"))
        run["manifest_checksums"]["output_manifest.tsv"] = hashlib.sha256(output_path.read_bytes()).hexdigest()
        run_path.write_text(json.dumps(run, sort_keys=True, indent=2) + "\n", encoding="utf-8")
        completed = subprocess.run([sys.executable, str(VALIDATE), str(self.report)], capture_output=True, text=True)
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("duplicate output identifier", completed.stderr)


if __name__ == "__main__":
    unittest.main()
