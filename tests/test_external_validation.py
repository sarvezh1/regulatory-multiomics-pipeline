import csv
import importlib.util
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


BENCHMARK = load_module("run_crispr_benchmark", ROOT / "bin" / "run_crispr_benchmark.py")


class ExternalValidationUnitTests(unittest.TestCase):
    def test_strict_overlap_excludes_bookended_interval(self):
        rows = [{"chromosome": "chr1", "start": "10", "end": "20", "region_id": "r1"}]
        index, maximum = BENCHMARK.interval_index(rows)
        self.assertEqual(BENCHMARK.overlaps(index, maximum, "chr1", 20, 30), [])
        self.assertEqual(BENCHMARK.overlaps(index, maximum, "chr1", 19, 30)[0]["region_id"], "r1")

    def test_binary_metrics_are_not_rank_metrics(self):
        records = [
            {"external_positive": True, "predicted": True},
            {"external_positive": True, "predicted": False},
            {"external_positive": False, "predicted": True},
            {"external_positive": False, "predicted": False},
        ]
        row = BENCHMARK.metric_row("predicted", 100000, "TEST", records)
        self.assertEqual((row["TP"], row["FP"], row["TN"], row["FN"]), ("1", "1", "1", "1"))
        self.assertEqual(row["MCC"], "0")

    def test_source_registry_declares_independence(self):
        path = ROOT / "metadata" / "external_validation_sources.tsv"
        with path.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        classes = {row["source_name"]: row["independence_class"] for row in rows}
        self.assertEqual(classes["ENCODE_K562_CCRE"], "ANNOTATION_CONCORDANCE_NONINDEPENDENT")
        self.assertEqual(classes["EPCrisprBenchmark_K562_TRAINING"], "FUNCTIONAL_INDEPENDENT")
        self.assertTrue(all(row["genome_build"].startswith("GRCh38") for row in rows))


if __name__ == "__main__":
    unittest.main()
