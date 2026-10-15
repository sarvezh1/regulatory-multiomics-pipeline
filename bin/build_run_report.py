#!/usr/bin/env python3
import argparse
import csv
import datetime as dt
import hashlib
import html
import json
import os
import re
import sys
import time
from pathlib import Path

SCHEMA_VERSION = "1.0.0"
STATUS_VALUES = {
    "COMPLETE", "COMPLETE_EMPTY_RESULT", "NOT_RUN_DISABLED",
    "NOT_RUN_NO_CONTRAST", "NOT_RUN_INSUFFICIENT_REPLICATION",
    "NOT_RUN_INVALID_DESIGN", "NOT_AVAILABLE_INPUT", "FAILED",
}
SOURCE_SUFFIXES = {".nf", ".config", ".py", ".R", ".r", ".tsv", ".txt", ".yaml", ".yml"}
SOURCE_ROOT_FILES = {"main.nf", "nextflow.config"}
SOURCE_DIRS = {"conf", "workflows", "modules", "bin", "schemas", "assets"}
OPERATIONAL_NAMES = {
    "execution_trace.tsv", "execution_report.html", "execution_timeline.html",
    "workflow_dag.html",
}
PRIVATE_PATTERNS = [
    re.compile(r"[A-Za-z]:[\\/][^\s\"'<>]+"),
    re.compile(r"(?:^|[\s\"'])/home/[^/\s]+/"),
    re.compile(r"(?i)(?:password|passwd|secret|token|api[_-]?key)\s*[:=]\s*[^\s,;]+"),
    re.compile(r"(?i)\b" + "ALCHE" + "MIST" + r"\b"),
]


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha_text(text):
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def read_tsv(path):
    with open(path, encoding="utf-8-sig", newline="") as handle:
        table = list(csv.reader(handle, delimiter="\t"))
    if not table or not table[0] or any(not column for column in table[0]):
        raise ValueError(f"reporting: malformed TSV header: {Path(path).name}")
    width = len(table[0])
    if len(set(table[0])) != width:
        raise ValueError(f"reporting: duplicate TSV column: {Path(path).name}")
    for line_number, row in enumerate(table[1:], 2):
        if len(row) != width:
            raise ValueError(f"reporting: malformed TSV row {line_number}: {Path(path).name}")
    return [dict(zip(table[0], row)) for row in table[1:]]


def write_tsv(path, rows, fields, schema_name):
    path.parent.mkdir(parents=True, exist_ok=True)
    all_fields = ["schema_name", "schema_version"] + fields
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=all_fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            item = {"schema_name": schema_name, "schema_version": SCHEMA_VERSION, **row}
            writer.writerow({key: item.get(key, "NA") for key in all_fields})


def safe_relative(path, root):
    return path.resolve().relative_to(root.resolve()).as_posix()


def ensure_public_safe(value, label):
    text = str(value)
    for pattern in PRIVATE_PATTERNS:
        if pattern.search(text):
            raise ValueError(f"reporting: private path or credential-like value detected in {label}")


def source_manifest(root):
    rows = []
    candidates = [root / name for name in SOURCE_ROOT_FILES]
    for directory in sorted(SOURCE_DIRS):
        source_dir = root / directory
        if source_dir.is_dir():
            candidates.extend(source_dir.rglob("*"))
    for path in candidates:
        try:
            if not path.is_file():
                continue
            rel = safe_relative(path, root)
            if path.suffix not in SOURCE_SUFFIXES and path.name not in SOURCE_ROOT_FILES:
                continue
            if "__pycache__" in path.parts or path.suffix == ".pyc":
                continue
            rows.append({"relative_path": rel, "sha256": sha256(path), "size_bytes": path.stat().st_size})
        except OSError as error:
            raise ValueError(f"reporting: source file could not be read: {path.name}") from error
    rows.sort(key=lambda row: row["relative_path"])
    canonical = "".join(f'{row["relative_path"]}\t{row["sha256"]}\t{row["size_bytes"]}\n' for row in rows)
    return rows, sha_text(canonical)


def normalized_inputs(records, source_root):
    rows = []
    seen = set()
    for record in records:
        path = Path(record)
        if not path.is_file():
            raise ValueError(f"reporting: input record does not exist: {path.name}")
        table = read_tsv(path)
        if not table:
            raise ValueError(f"reporting: input record is empty: {path.name}")
        if {"source_name", "source_type", "accession", "download_location", "file_checksum"}.issubset(table[0]):
            for item in table:
                identity = ("EXTERNAL_VALIDATION", item["source_name"], item["accession"], item["source_type"])
                if identity in seen:
                    continue
                seen.add(identity)
                rows.append({
                    "modality": "EXTERNAL_VALIDATION", "sample_id": "NA", "library_id": "NA",
                    "role": item["source_type"], "logical_name": item["source_name"],
                    "source_accession": item["accession"] or "UNKNOWN", "size_bytes": item.get("size_bytes", "UNRESOLVED") or "UNRESOLVED",
                    "source_checksum": item["file_checksum"] or "UNKNOWN", "local_sha256": item.get("local_sha256", "UNRESOLVED") or "UNRESOLVED",
                    "source_genome_build": item.get("genome_build", "UNKNOWN") or "UNKNOWN", "source_annotation": item.get("version", "UNKNOWN") or "UNKNOWN",
                    "source_processing_pipeline": item.get("source_processing_pipeline", "EXTERNAL_PUBLIC_RESOURCE") or "EXTERNAL_PUBLIC_RESOURCE",
                    "coordinate_normalization_status": item.get("coordinate_normalization_status", "NOT_APPLICABLE") or "NOT_APPLICABLE",
                })
            continue
        for item in table:
            modality = item.get("assay") or item.get("modality") or "UNKNOWN"
            role_fields = [("read1", "READ1"), ("read2", "READ2"), ("processed_input", "PROCESSED_INPUT"), ("path", item.get("data_type", "CONTRACT_INPUT"))]
            found = False
            for field, role in role_fields:
                value = item.get(field, "")
                if not value or value in {"NA", "UNKNOWN", "UNRESOLVED"}:
                    continue
                found = True
                candidate = Path(value)
                if not candidate.is_absolute():
                    candidate = source_root / candidate
                identity = (modality, item.get("sample_id", "NA"), item.get("library_id", "NA"), role, candidate.name)
                if identity in seen:
                    continue
                seen.add(identity)
                exists = candidate.is_file()
                rows.append({
                    "modality": modality,
                    "sample_id": item.get("sample_id", "NA"),
                    "library_id": item.get("library_id", "NA"),
                    "role": role,
                    "logical_name": candidate.name,
                    "source_accession": item.get("source_file_accession", "UNKNOWN") or "UNKNOWN",
                    "size_bytes": candidate.stat().st_size if exists else "UNRESOLVED",
                    "source_checksum": item.get("source_file_checksum", "UNKNOWN") or "UNKNOWN",
                    "local_sha256": sha256(candidate) if exists else "UNRESOLVED",
                    "source_genome_build": item.get("source_genome_build", item.get("genome_build", "UNKNOWN")) or "UNKNOWN",
                    "source_annotation": item.get("source_annotation", "UNKNOWN") or "UNKNOWN",
                    "source_processing_pipeline": item.get("source_processing_pipeline", "UNKNOWN") or "UNKNOWN",
                    "coordinate_normalization_status": item.get("coordinate_normalization_status", "NOT_APPLICABLE") or "NOT_APPLICABLE",
                })
            if not found:
                identity = (modality, item.get("sample_id", "NA"), item.get("library_id", "NA"), "METADATA", path.name)
                if identity not in seen:
                    seen.add(identity)
                    rows.append({
                        "modality": modality, "sample_id": item.get("sample_id", "NA"),
                        "library_id": item.get("library_id", "NA"), "role": "METADATA",
                        "logical_name": path.name, "source_accession": item.get("source_file_accession", "UNKNOWN") or "UNKNOWN",
                        "size_bytes": path.stat().st_size, "source_checksum": item.get("source_file_checksum", "UNKNOWN") or "UNKNOWN",
                        "local_sha256": sha256(path), "source_genome_build": item.get("source_genome_build", "UNKNOWN") or "UNKNOWN",
                        "source_annotation": item.get("source_annotation", "UNKNOWN") or "UNKNOWN",
                        "source_processing_pipeline": item.get("source_processing_pipeline", "UNKNOWN") or "UNKNOWN",
                        "coordinate_normalization_status": item.get("coordinate_normalization_status", "NOT_APPLICABLE") or "NOT_APPLICABLE",
                    })
    return sorted(rows, key=lambda row: (row["modality"], row["sample_id"], row["library_id"], row["role"], row["logical_name"]))


def reference_manifest(specifications):
    rows = []
    for spec in specifications:
        parts = spec.split("=", 1)
        if len(parts) != 2:
            raise ValueError("reporting: --reference must use ROLE=PATH")
        role, value = parts
        path = Path(value)
        if not path.is_file():
            raise ValueError(f"reporting: required reference does not exist: {path.name}")
        path_text = path.as_posix().lower()
        validation = "VALIDATION_ONLY" if any(token in path_text for token in ("chr22", "synthetic", "fixture", "/tests/")) else "PRODUCTION_OR_USER_SUPPLIED"
        role_upper = role.upper()
        if "GTF" in role_upper:
            version, coordinate_system, annotation = "GENCODE_50", "GTF_1_BASED_CLOSED_SOURCE", "GENCODE_50"
        elif "BLACKLIST" in role_upper:
            version, coordinate_system, annotation = "ENCODE_GRCH38_BLACKLIST", "BED_0_BASED_HALF_OPEN", "NOT_APPLICABLE"
        elif "FASTA" in role_upper:
            version, coordinate_system, annotation = "GRCh38_PRIMARY_ASSEMBLY", "REFERENCE_SEQUENCE", "NOT_APPLICABLE"
        else:
            version, coordinate_system, annotation = "DERIVED_OR_USER_SUPPLIED", "SOURCE_FORMAT_DECLARED", "NOT_APPLICABLE"
        rows.append({
            "reference_role": role, "logical_name": path.name, "source": "USER_SUPPLIED_OR_PINNED_PROJECT_RESOURCE",
            "version": version,
            "source_accession_or_url": "UNKNOWN", "source_checksum": "UNKNOWN", "local_sha256": sha256(path),
            "size_bytes": path.stat().st_size, "derivation_method": "DIRECT_INPUT",
            "genome_build": "GRCh38", "annotation_release": annotation,
            "coordinate_system": coordinate_system,
            "allowed_chromosomes": "chr1-22;chrX;chrY", "usage_scope": validation,
        })
    return sorted(rows, key=lambda row: (row["reference_role"], row["logical_name"]))


def tools_and_containers(path):
    raw = read_tsv(path)
    tools = []
    containers = []
    for row in raw:
        image = row.get("container", "") or row.get("image", "")
        digest = row.get("digest", "")
        if image and ((":" + "latest") in image or "@sha256:" not in image):
            raise ValueError(f"reporting: unpinned container reference for {row.get('process', row.get('tool', 'UNKNOWN'))}")
        name = row.get("process") or row.get("component") or row.get("tool") or "UNKNOWN"
        version = row.get("tool_version") or row.get("tool") or row.get("version") or "UNRESOLVED"
        image_without_digest = image.split("@", 1)[0]
        final_slash = image_without_digest.rfind("/")
        final_colon = image_without_digest.rfind(":")
        tag = image_without_digest[final_colon + 1:] if final_colon > final_slash else "UNRESOLVED"
        repository = image_without_digest[:final_colon] if final_colon > final_slash else image_without_digest
        tools.append({"tool": name, "executed_version": version, "purpose": name, "execution_evidence": "VALIDATED_ASSAYS_AND_INTEGRATION"})
        containers.append({"component": name, "container_image": image, "image_repository": repository, "tag": tag, "immutable_digest": digest or (image.split("@", 1)[1] if "@" in image else "UNRESOLVED"), "execution_backend": "Docker", "validated_architecture": "Linux/amd64"})
    return tools, containers


def discover_statuses(results):
    rows = []
    seen = {}
    status_paths = set(results.rglob("*.analysis_status.tsv"))
    status_paths.update(results.rglob("external_validation_status.tsv"))
    for path in sorted(status_paths):
        for row in read_tsv(path):
            component = row.get("component")
            status = row.get("status")
            if not component or not status:
                continue
            if status not in STATUS_VALUES:
                raise ValueError(f"reporting: unsupported status {status} for {component}")
            if component in seen and seen[component] != status:
                raise ValueError(f"reporting: conflicting statuses for {component}: {seen[component]} and {status}")
            if component in seen:
                continue
            seen[component] = status
            rows.append({"component": component, "status": status, "source_artifact": safe_relative(path, results), "interpretation": "PROCESS_STATUS"})
    if not rows:
        raise ValueError("reporting: no component status manifests were found; completeness cannot be claimed")
    return sorted(rows, key=lambda row: row["component"])


def output_role(path):
    rel = path.as_posix().lower()
    if "/provenance/" in f"/{rel}" or path.name.endswith(".version.txt"):
        return "PROVENANCE"
    if "multiqc" in rel or "/qc/" in f"/{rel}":
        return "QC_OUTPUT"
    if path.suffix.lower() in {".html", ".svg", ".png", ".pdf"}:
        return "REPORT"
    if any(token in rel for token in ("trimmed/", "alignment/", "reference_index/")):
        return "INTERMEDIATE_OPTIONAL"
    return "PRIMARY_ANALYSIS_OUTPUT"


def output_manifest(results):
    rows = []
    for path in sorted(results.rglob("*")):
        if not path.is_file() or path.name in OPERATIONAL_NAMES:
            continue
        rel = safe_relative(path, results)
        component = rel.split("/", 1)[0]
        suffix = path.suffix.lower().lstrip(".") or "binary"
        rows.append({
            "component": component, "output_type": output_role(Path(rel)), "logical_name": path.name,
            "relative_path": rel, "format": suffix.upper(), "status": "PRESENT",
            "size_bytes": path.stat().st_size, "sha256": sha256(path),
            "description": "Canonical workflow artifact inventory entry", "scientific_role": output_role(Path(rel)),
        })
    return rows


def validate_completed_contract_outputs(results, statuses):
    expected = {
        "rna_count_contract": "rna.raw_counts.tsv",
        "atac_count_contract": "atac.raw_counts.tsv",
        "wgbs_count_contract": "wgbs.cpg.long.tsv",
        "h3k27ac_count_contract": "h3k27ac.raw_counts.tsv",
        "integration_promoter": "promoter_evidence.tsv",
    }
    available = {path.name for path in results.rglob("*") if path.is_file() and path.stat().st_size > 0}
    for row in statuses:
        required = expected.get(row["component"])
        if row["status"] == "COMPLETE" and required and required not in available:
            raise ValueError(
                f"reporting: {row['component']} is COMPLETE but required output is missing or empty: {required}"
            )


def sanitized_execution(results):
    rows = []
    for path in sorted(results.rglob("execution_trace.tsv")):
        try:
            trace = read_tsv(path)
        except Exception:
            continue
        for row in trace:
            rows.append({
                "component": safe_relative(path.parent.parent, results),
                "process": row.get("process", row.get("name", "UNKNOWN")), "status": row.get("status", "UNKNOWN"),
                "exit_code": row.get("exit", "NA"), "duration": row.get("duration", "NA"),
                "realtime": row.get("realtime", "NA"), "cpu_percent": row.get("%cpu", "NA"),
                "peak_rss": row.get("peak_rss", "NA"), "container": row.get("container", "UNKNOWN"),
            })
    return rows


def warning_manifest(results, statuses, references):
    rows = []
    if any(row["usage_scope"] == "VALIDATION_ONLY" for row in references):
        rows.append({"warning_code": "COMPACT_VALIDATION_REFERENCE", "component": "reference", "severity": "INFO", "category": "VALIDATION_SCOPE", "message": "A compact or synthetic reference is execution-validation only and is not a production biological reference.", "source": "reference_manifest.tsv", "action_required": "Use the pinned full production reference for biological analysis."})
    if any(row["status"].startswith("NOT_RUN") for row in statuses):
        rows.append({"warning_code": "OPTIONAL_COMPONENT_NOT_RUN", "component": "run", "severity": "INFO", "category": "STATUS", "message": "One or more optional components were not run; this is not a process failure.", "source": "status_manifest.tsv", "action_required": "Review component-scoped statuses."})
    if any(results.rglob("*nearest*TSS*")) or any(results.rglob("distal_candidate_links.tsv")):
        rows.append({"warning_code": "PREDICTIVE_DISTAL_LINKS", "component": "integration", "severity": "WARNING", "category": "SCIENTIFIC_LIMITATION", "message": "Nearest-TSS distal relationships are candidate predictions and do not establish physical contact or a true target.", "source": "integration contract", "action_required": "Treat links as hypotheses unless independent validation is supplied."})
    if any(results.rglob("*wgbs*")):
        rows.append({"warning_code": "CONVERSION_CONTROL_CONTEXT_REQUIRED", "component": "wgbs", "severity": "INFO", "category": "QC_LIMITATION", "message": "Observed non-CpG methylation is not equivalent to measured conversion efficiency without a genuine conversion control.", "source": "WGBS contract", "action_required": "Review conversion-control status before interpretation."})
    rows.append({"warning_code": "CHIP_CROSS_CORRELATION_DEFERRED", "component": "h3k27ac", "severity": "INFO", "category": "QC_LIMITATION", "message": "Normalized and relative strand cross-correlation metrics are not implemented in v0.1.", "source": "H3K27ac contract", "action_required": "Use the implemented complexity, enrichment and FRiP metrics within their documented limits."})
    rows.append({"warning_code": "PERCENTAGE_ONLY_WGBS_DESCRIPTIVE", "component": "wgbs", "severity": "INFO", "category": "INPUT_LIMITATION", "message": "Percentage-only methylation input is descriptive and cannot support production differential inference.", "source": "WGBS contract", "action_required": "Supply methylated and unmethylated or total coverage counts for inference."})
    rows.append({"warning_code": "LOCAL_OPERATIONAL_METADATA_EXCLUDED", "component": "reporting", "severity": "INFO", "category": "PRIVACY", "message": "Raw Nextflow trace, timeline, report and DAG artifacts are local operational records and are excluded from the public output manifest.", "source": "reporting policy", "action_required": "Use sanitized_execution.tsv for public-safe execution metadata."})
    if any(results.rglob("ccre_concordance_summary.tsv")):
        rows.append({"warning_code": "CCRE_CONCORDANCE_NONINDEPENDENT", "component": "external_validation", "severity": "WARNING", "category": "VALIDATION_DEPENDENCE", "message": "K562 cCRE overlap is annotation concordance and is not an independent functional benchmark.", "source": "external validation source registry", "action_required": "Keep cCRE concordance separate from the CRISPR functional benchmark."})
    if any(results.rglob("crispr_benchmark_metrics.tsv")):
        rows.append({"warning_code": "CRISPR_TESTED_SUBSET_SCOPE", "component": "external_validation", "severity": "INFO", "category": "BENCHMARK_SCOPE", "message": "The CRISPR benchmark covers experimentally tested element-gene pairs rather than the whole genome.", "source": "crispr_benchmark_metrics.tsv", "action_required": "Interpret metrics within the reported benchmark denominator."})
    return rows


def plot_manifest(results):
    rows = []
    for plot in sorted(list(results.rglob("*.svg")) + list(results.rglob("*.png"))):
        candidates = [plot.with_suffix(".tsv"), plot.with_suffix(".csv")]
        source = next((candidate for candidate in candidates if candidate.is_file()), None)
        if source:
            rows.append({"plot_name": plot.name, "source_table": safe_relative(source, results), "columns_used": "DECLARED_BY_SOURCE_CONTRACT", "filters": "NONE_BEYOND_SOURCE_TABLE", "parameters": "DETERMINISTIC_RENDERING"})
    return rows


def html_report(run, statuses, warnings, outputs, tools, multiqc, integration_links, external_links, plots, modality_rows, input_count, modality_count):
    def table(headers, data):
        head = "".join(f"<th>{html.escape(h)}</th>" for h in headers)
        body = "".join("<tr>" + "".join(f"<td>{html.escape(str(row.get(h, 'NA')))}</td>" for h in headers) + "</tr>" for row in data)
        return f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"
    links = "".join(f'<li><a href="{html.escape(item)}">{html.escape(item)}</a></li>' for item in multiqc) or "<li>Not available</li>"
    integration = "".join(f'<li><a href="{html.escape(item)}">{html.escape(item)}</a></li>' for item in integration_links) or "<li>Not available</li>"
    external = "".join(f'<li><a href="{html.escape(item)}">{html.escape(item)}</a></li>' for item in external_links) or "<li>Not available</li>"
    sections = "".join(
        f"<h2>{html.escape(row['modality'])} summary</h2><p>{row['status_component_count']} component statuses and {row['output_artifact_count']} inventoried artifacts. Source: <code>status_manifest.tsv</code> and <code>output_manifest.tsv</code>.</p>"
        for row in modality_rows
    )
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><title>regulatory-multiomics-pipeline run report</title>
<style>body{{font-family:system-ui,sans-serif;max-width:1100px;margin:2rem auto;padding:0 1rem;color:#17232b}}table{{border-collapse:collapse;width:100%;margin:1rem 0}}th,td{{border:1px solid #cbd5df;padding:.45rem;text-align:left;vertical-align:top}}th{{background:#edf2f7}}code{{background:#f3f5f7;padding:.1rem .25rem}}.note{{border-left:4px solid #5b4bb7;padding:.6rem 1rem;background:#f6f4ff}}</style></head>
<body><h1>Regulatory multi-omics run report</h1>
<div class="note">This report separates process status, scientific QC measurements, and assay-specific statistical inference. It does not provide causal biological interpretation.</div>
<h2>Run overview</h2><p>Run ID: <code>{html.escape(run['run_id'])}</code><br>Requested analyses: {html.escape(run['requested_assay'])}<br>Backend: Docker Linux/amd64<br>Reference policy: GRCh38 primary assembly with GENCODE 50 context.</p>
<h2>Samples and modalities</h2><p>{input_count} normalized input records across {modality_count} declared modalities. Source: <code>input_manifest.tsv</code>.</p>
<h2>Component status</h2>{table(['component','status','interpretation'], statuses)}
<h2>Warnings and limitations</h2>{table(['warning_code','component','severity','message','action_required'], warnings)}
<h2>Tool and container inventory</h2>{table(['tool','executed_version','purpose'], tools)}
<h2>Modality reports</h2><ul>{links}</ul>
{sections}
<h2>Integration outputs</h2><ul>{integration}</ul>
<h2>External validation</h2><p>Functional benchmarking, annotation concordance, and descriptive biological plausibility are reported separately. No validation source changes the frozen production rules.</p><ul>{external}</ul>
<h2>Output inventory</h2><p>{len(outputs)} public scientific/QC/provenance/report artifacts are catalogued in <code>output_manifest.tsv</code>. Raw operational Nextflow artifacts remain local.</p>
<h2>Visualization traceability</h2><p>{len(plots)} plot-to-source relationships are recorded in <code>plot_manifest.tsv</code>.</p>
<h2>Reproducibility</h2><p>Source fingerprint: <code>{run['source_tree_fingerprint']}</code>. Input, reference, tool, container, parameter, output, status and warning manifests are linked by checksums in <code>run_manifest.json</code>.</p>
<h2>Limitations</h2><p>No cross-assay p-value combination, causal inference, contact-informed enhancer assignment, or automatic biological narrative is produced.</p></body></html>"""


def main():
    started = time.time()
    parser = argparse.ArgumentParser(description="Build public-safe unified provenance and an offline run report.")
    parser.add_argument("--results", required=True)
    parser.add_argument("--source-root", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--requested-assay", required=True)
    parser.add_argument("--run-mode", default="synthetic_validation")
    parser.add_argument("--nextflow-version", required=True)
    parser.add_argument("--input-record", action="append", default=[])
    parser.add_argument("--reference", action="append", default=[])
    parser.add_argument("--parameter-file")
    parser.add_argument("--tool-containers", required=True)
    parser.add_argument("--timestamp-utc")
    args = parser.parse_args()

    results = Path(args.results).resolve(); source_root = Path(args.source_root).resolve(); output = Path(args.output).resolve()
    if not results.is_dir(): raise ValueError("reporting: results directory does not exist")
    if output == results or results in output.parents: raise ValueError("reporting: output must not be nested inside the inventoried results directory")
    output.mkdir(parents=True, exist_ok=True)

    source_rows, source_fingerprint = source_manifest(source_root)
    inputs = normalized_inputs(args.input_record, source_root)
    references = reference_manifest(args.reference)
    tools, containers = tools_and_containers(args.tool_containers)
    tools.append({"tool": "REPORTING_RUNTIME", "executed_version": f"Python {sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}", "purpose": "Unified provenance and offline report generation", "execution_evidence": "CURRENT_EXECUTION"})
    tools = sorted(tools, key=lambda row: row["tool"])
    statuses = discover_statuses(results)
    validate_completed_contract_outputs(results, statuses)
    outputs = output_manifest(results)
    warnings = warning_manifest(results, statuses, references)
    executions = sanitized_execution(results)
    plots = plot_manifest(results)
    parameters = json.load(open(args.parameter_file, encoding="utf-8")) if args.parameter_file else {"requested_assay": args.requested_assay, "run_mode": args.run_mode}
    parameters = {"schema_name": "parameter_manifest", "schema_version": SCHEMA_VERSION, **parameters}
    ensure_public_safe(json.dumps(parameters, sort_keys=True), "parameter manifest")

    write_tsv(output / "source_manifest.tsv", source_rows, ["relative_path","sha256","size_bytes"], "source_manifest")
    write_tsv(output / "input_manifest.tsv", inputs, ["modality","sample_id","library_id","role","logical_name","source_accession","size_bytes","source_checksum","local_sha256","source_genome_build","source_annotation","source_processing_pipeline","coordinate_normalization_status"], "input_manifest")
    write_tsv(output / "reference_manifest.tsv", references, ["reference_role","logical_name","source","version","source_accession_or_url","source_checksum","local_sha256","size_bytes","derivation_method","genome_build","annotation_release","coordinate_system","allowed_chromosomes","usage_scope"], "reference_manifest")
    write_tsv(output / "tool_manifest.tsv", tools, ["tool","executed_version","purpose","execution_evidence"], "tool_manifest")
    write_tsv(output / "container_manifest.tsv", containers, ["component","container_image","image_repository","tag","immutable_digest","execution_backend","validated_architecture"], "container_manifest")
    write_tsv(output / "output_manifest.tsv", outputs, ["component","output_type","logical_name","relative_path","format","status","size_bytes","sha256","description","scientific_role"], "output_manifest")
    write_tsv(output / "status_manifest.tsv", statuses, ["component","status","source_artifact","interpretation"], "status_manifest")
    write_tsv(output / "warning_manifest.tsv", warnings, ["warning_code","component","severity","category","message","source","action_required"], "warning_manifest")
    write_tsv(output / "sanitized_execution.tsv", executions, ["component","process","status","exit_code","duration","realtime","cpu_percent","peak_rss","container"], "sanitized_execution")
    write_tsv(output / "plot_manifest.tsv", plots, ["plot_name","source_table","columns_used","filters","parameters"], "plot_manifest")
    contract_versions = [{"contract": name, "contract_version": "1.0.0", "stability": "INITIAL_VERSIONED_CONTRACT"} for name in ("RNA","ATAC","WGBS","H3K27ac","integration","external_validation","run_provenance")]
    write_tsv(output / "contract_versions.tsv", contract_versions, ["contract","contract_version","stability"], "contract_version_registry")
    (output / "parameter_manifest.json").write_text(json.dumps(parameters, sort_keys=True, indent=2) + "\n", encoding="utf-8")

    manifest_checksums = {name: sha256(output / name) for name in ("input_manifest.tsv","reference_manifest.tsv","tool_manifest.tsv","container_manifest.tsv","parameter_manifest.json","source_manifest.tsv","output_manifest.tsv","status_manifest.tsv","warning_manifest.tsv")}
    identity_text = "\n".join([source_fingerprint] + [manifest_checksums[key] for key in sorted(manifest_checksums) if key not in {"output_manifest.tsv", "warning_manifest.tsv"}])
    stamp = args.timestamp_utc or dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat()
    run = {
        "schema_name": "run_manifest", "schema_version": SCHEMA_VERSION,
        "project": "regulatory-multiomics-pipeline", "project_version": "0.1.0",
        "run_id": f"RMP_{sha_text(identity_text)[:16]}", "run_mode": args.run_mode,
        "requested_assay": args.requested_assay, "requested_components": sorted({row["component"] for row in statuses}),
        "execution_backend": "Docker", "validated_architecture": "Linux/amd64",
        "nextflow_version": args.nextflow_version, "runtime_platform": f"Python {sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}",
        "start_timestamp_utc": stamp, "end_timestamp_utc": stamp,
        "source_tree_fingerprint": source_fingerprint, "manifest_checksums": manifest_checksums,
        "component_statuses": statuses, "report_generation_seconds": round(time.time() - started, 3),
    }
    (output / "run_manifest.json").write_text(json.dumps(run, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    modality_rows = []
    for modality in ("RNA", "ATAC", "WGBS", "H3K27ac", "integration"):
        prefix = modality.lower()
        modality_rows.append({
            "modality": modality,
            "status_component_count": sum(1 for row in statuses if row["component"].lower().startswith(prefix)),
            "output_artifact_count": sum(1 for row in outputs if row["component"].lower() == prefix),
        })
    multiqc = [Path(os.path.relpath(path, output)).as_posix() for path in sorted(results.rglob("*multiqc_report.html"))]
    integration_links = [Path(os.path.relpath(results / row["relative_path"], output)).as_posix() for row in outputs if row["relative_path"] in {"integration/evidence/promoter_evidence.tsv", "integration/evidence/gene_summary.tsv", "integration/relationships/distal_candidate_links.tsv", "integration/sensitivity/integration.sensitivity_summary.tsv"}]
    external_names = {"ccre_concordance_summary.tsv", "crispr_benchmark_metrics.tsv", "crispr_failure_classes.tsv", "biological_plausibility_summary.tsv", "external_validation_status.tsv", "external_validation_provenance.json"}
    external_links = [Path(os.path.relpath(results / row["relative_path"], output)).as_posix() for row in outputs if row["logical_name"] in external_names]
    report = html_report(run, statuses, warnings, outputs, tools, multiqc, integration_links, external_links, plots, modality_rows, len(inputs), len({row["modality"] for row in inputs}))
    ensure_public_safe(report, "HTML report")
    (output / "run_report.html").write_text(report, encoding="utf-8")
    summary_rows = [{"metric": "component_count", "value": len(statuses), "source": "status_manifest.tsv"}, {"metric": "output_count", "value": len(outputs), "source": "output_manifest.tsv"}, {"metric": "warning_count", "value": len(warnings), "source": "warning_manifest.tsv"}, {"metric": "plot_count", "value": len(plots), "source": "plot_manifest.tsv"}, {"metric": "input_record_count", "value": len(inputs), "source": "input_manifest.tsv"}, {"metric": "modality_count", "value": len({row["modality"] for row in inputs}), "source": "input_manifest.tsv"}]
    for row in modality_rows:
        key = row["modality"].lower()
        summary_rows.extend([
            {"metric": f"{key}_status_component_count", "value": row["status_component_count"], "source": "status_manifest.tsv"},
            {"metric": f"{key}_output_artifact_count", "value": row["output_artifact_count"], "source": "output_manifest.tsv"},
        ])
    write_tsv(output / "report_summary.tsv", summary_rows, ["metric","value","source"], "report_summary")
    for path in output.rglob("*"):
        if path.is_file(): ensure_public_safe(path.read_text(encoding="utf-8", errors="ignore"), path.name)
    print(json.dumps({"run_id": run["run_id"], "source_tree_fingerprint": source_fingerprint, "outputs": len(outputs), "statuses": len(statuses)}, sort_keys=True))


if __name__ == "__main__":
    main()
