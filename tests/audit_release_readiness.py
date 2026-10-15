#!/usr/bin/env python3
import argparse
import csv
import json
import os
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
VERSION = "0.1.0"
AUTHOR = "Sarvesh Galgale"

parser = argparse.ArgumentParser(description="Audit future-public project files before release.")
parser.add_argument("--allow-git", action="store_true", help="Allow the checkout metadata present in hosted CI.")
args = parser.parse_args()

REQUIRED = [
    "README.md", "LICENSE", "CITATION.cff", "main.nf", "nextflow.config", ".gitignore", ".gitattributes",
    ".github/workflows/ci.yml", "docs/documentation_map.md", "docs/public_file_inventory.md",
    "docs/reporting.md", "docs/reproducibility.md", "docs/output_reference.md",
    "docs/usage.md", "docs/public_data_sources.md", "docs/external_validation.md", "CHANGELOG.md",
    "docs/validation/rna.md", "docs/validation/atac.md", "docs/validation/wgbs.md",
    "docs/validation/h3k27ac.md", "docs/validation/integration.md",
    "docs/validation/reproducibility.md", "docs/validation/external.md",
    "assets/readme/external_benchmark_metrics.svg", "assets/readme/false_negative_breakdown.svg",
    "assets/readme/k562_promoter_evidence.svg",
    "metadata/external_validation_sources.tsv", "metadata/external_validation_prediction_freeze.tsv",
    "schemas/run_provenance_contract.tsv", "schemas/contract_versions.tsv",
]
PUBLIC_ROOTS = [".github", "bin", "conf", "workflows", "modules", "schemas", "assets", "docs", "metadata", "tests"]
TEXT_SUFFIXES = {".cff", ".md", ".nf", ".config", ".py", ".R", ".r", ".tsv", ".txt", ".yaml", ".yml", ".ps1", ".json"}
EXCLUDED_PARTS = {"output", "work", "__pycache__", "generated"}
PROHIBITED = [
    "Cod" + "ex", "Chat" + "GPT", "Open" + "AI", "AI-" + "assisted",
    "AI " + "assistance", "generated " + "by", "assist" + "ant", "assist" + "ants",
    "pro" + "mpt", "pro" + "mpts", "ag" + "ent", "ag" + "ents", "language " + "model",
]
PRIVATE = [
    re.compile(r"[A-Za-z]:[\\/](?:Users|Pipelines)[\\/]", re.I),
    re.compile(r"(?:^|[\s\"'])/home/[^/\s]+/"),
    re.compile(r"\b" + "ALCHE" + "MIST" + r"\b", re.I),
    re.compile(r"(?:password|secret|token|api[_-]?key)\s*[:=]\s*[\"'][^\"']{8,}[\"']", re.I),
]

issues = []
for relative in REQUIRED:
    path = ROOT / relative
    if not path.is_file() or path.stat().st_size == 0:
        issues.append(f"missing required file: {relative}")
if (ROOT / ".git").exists() and not args.allow_git:
    issues.append("a Git repository exists before the authorized initialization phase")

files = [ROOT / name for name in ("README.md", "LICENSE", "CITATION.cff", "main.nf", "nextflow.config", ".gitignore", ".gitattributes")]
for directory in PUBLIC_ROOTS:
    search_root = ROOT / directory
    if not search_root.exists():
        continue
    for current, dirs, names in os.walk(search_root, topdown=True):
        dirs[:] = [name for name in dirs if name not in EXCLUDED_PARTS]
        if Path(current).resolve() == (ROOT / "tests" / "integration").resolve():
            dirs[:] = []
        for name in names:
            path = Path(current) / name
            if path.suffix not in TEXT_SUFFIXES:
                continue
            if path.name.endswith(".validated.tsv"):
                continue
            files.append(path)

for path in sorted(set(files)):
    text = path.read_text(encoding="utf-8", errors="ignore")
    for term in PROHIBITED:
        if re.search(rf"\b{re.escape(term)}\b", text, re.I):
            issues.append(f"{path.relative_to(ROOT)}: prohibited public-project term")
    for pattern in PRIVATE:
        if pattern.search(text):
            issues.append(f"{path.relative_to(ROOT)}: private path, machine identifier, or credential-like value")
    if (":" + "latest") in text:
        issues.append(f"{path.relative_to(ROOT)}: floating container tag")
    unresolved_terms = ["YOUR " + "NAME", "T" + "BD", "maintainer " + "placeholder", "author " + "placeholder"]
    if any(re.search(rf"\b{re.escape(term)}\b", text, re.I) for term in unresolved_terms):
        issues.append(f"{path.relative_to(ROOT)}: unresolved public metadata placeholder")

markdown = [path for path in files if path.suffix == ".md"]
link_pattern = re.compile(r"\[[^\]]+\]\(([^)]+)\)")
for path in markdown:
    text = path.read_text(encoding="utf-8", errors="ignore")
    for target in link_pattern.findall(text):
        if target.startswith(("http://", "https://", "#", "mailto:")):
            continue
        destination = (path.parent / target.split("#", 1)[0]).resolve()
        if not destination.exists():
            issues.append(f"{path.relative_to(ROOT)}: broken local link: {target}")

readme = (ROOT / "README.md").read_text(encoding="utf-8")
config = (ROOT / "nextflow.config").read_text(encoding="utf-8")
citation = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
license_text = (ROOT / "LICENSE").read_text(encoding="utf-8")
workflow = (ROOT / ".github/workflows/ci.yml").read_text(encoding="utf-8")

for label, text in (("README.md", readme), ("nextflow.config", config), ("CITATION.cff", citation)):
    if VERSION not in text:
        issues.append(f"{label}: expected pre-release version {VERSION}")
if "Copyright (c) 2026 Sarvesh Galgale" not in license_text or not license_text.startswith("MIT License"):
    issues.append("LICENSE: standard MIT heading or copyright is incorrect")
for fragment in ("cff-version: 1.2.0", "title: regulatory-multiomics-pipeline", "type: software", "license: MIT", "family-names: Galgale", "given-names: Sarvesh"):
    if fragment not in citation:
        issues.append(f"CITATION.cff: missing required metadata: {fragment}")
if "flowchart TD" not in readme or "External validation" not in readme:
    issues.append("README.md: final vertical workflow or external-validation node is missing")
for path in (ROOT / "bin").glob("*"):
    if path.is_file() and re.search(r"0\.[0-9]+\.0-block[0-9]+", path.read_text(encoding="utf-8", errors="ignore")):
        issues.append(f"{path.relative_to(ROOT)}: active runtime metadata still uses a development-state version")
for fragment in ("actions/checkout@v4", "actions/setup-python@v5", "audit_release_readiness.py --allow-git", "create_synthetic_fixtures.py", "create_atac_fixtures.py", "create_wgbs_fixtures.py", "create_chip_fixtures.py", "create_integration_fixtures.py", "test_*contracts.py"):
    if fragment not in workflow:
        issues.append(f".github/workflows/ci.yml: missing required CI step: {fragment}")

ignore_text = (ROOT / ".gitignore").read_text(encoding="utf-8")
for fragment in ("work/", "results/", ".nextflow*", "*.fastq.gz", "*.bam", "*.bai", "*.bigWig", "*.env", "tests/integration/", "external_data/"):
    if fragment not in ignore_text:
        issues.append(f".gitignore: missing safety rule: {fragment}")

for path in sorted(set(files)):
    if path.suffix == ".json":
        try:
            json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError as exc:
            issues.append(f"{path.relative_to(ROOT)}: malformed JSON: {exc.msg}")
    if path.suffix == ".tsv":
        try:
            with path.open(encoding="utf-8", newline="") as handle:
                rows = [row for row in csv.reader(handle, delimiter="\t") if any(field.strip() for field in row)]
            if not rows or not rows[0] or any(len(row) != len(rows[0]) for row in rows[1:]):
                issues.append(f"{path.relative_to(ROOT)}: malformed TSV column structure")
        except (OSError, UnicodeError) as exc:
            issues.append(f"{path.relative_to(ROOT)}: unreadable TSV: {exc}")

large_candidates = []
for relative in ("bin", "conf", "workflows", "modules", "schemas", "assets", "docs", "metadata"):
    for path in (ROOT / relative).rglob("*"):
        if path.is_file() and path.stat().st_size > 10 * 1024 * 1024:
            large_candidates.append(str(path.relative_to(ROOT)))
if large_candidates:
    issues.append("unexpected public-source file larger than 10 MiB: " + ", ".join(large_candidates))

if issues:
    raise SystemExit("\n".join(sorted(set(f"release audit: {item}" for item in issues))))
print(f"Release-readiness audit passed: {len(set(files))} public text artifacts checked; metadata, links, CI, hygiene, and ignore rules are coherent")
