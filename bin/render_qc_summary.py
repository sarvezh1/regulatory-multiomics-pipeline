#!/usr/bin/env python3
import argparse
import csv
import html
import shutil
from pathlib import Path


SAFE_PREFIXES = (
    "ataqv_",
    "bowtie2_",
    "fastp-",
    "fastp_",
    "fastqc-",
    "fastqc_",
    "multiqc_ataqv",
    "multiqc_bowtie2",
    "multiqc_fastp",
    "multiqc_fastqc",
    "multiqc_featurecounts",
    "multiqc_general_stats",
    "multiqc_samtools",
    "multiqc_software_versions",
    "multiqc_star",
    "samtools-",
    "samtools_",
    "star_",
)


def read_table(path: Path) -> tuple[list[str], list[list[str]]]:
    if not path.exists() or path.stat().st_size == 0:
        return [], []
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if not rows:
        return [], []
    return rows[0], rows[1:]


def table_markup(columns: list[str], rows: list[list[str]]) -> str:
    if not columns:
        return "<p>No aggregate metrics were available for this input route.</p>"
    head = "".join(f"<th>{html.escape(value)}</th>" for value in columns)
    body = []
    for row in rows:
        padded = row + [""] * max(0, len(columns) - len(row))
        cells = "".join(f"<td>{html.escape(value)}</td>" for value in padded[: len(columns)])
        body.append(f"<tr>{cells}</tr>")
    return f"<table><thead><tr>{head}</tr></thead><tbody>{''.join(body)}</tbody></table>"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--title", required=True)
    args = parser.parse_args()

    source_dir = args.input_dir
    if not (source_dir / "multiqc_general_stats.txt").exists():
        candidates = sorted(
            path
            for path in source_dir.iterdir()
            if path.is_dir()
            and (path.name == "multiqc_data" or path.name.endswith("_multiqc_report_data"))
        )
        source_dir = candidates[0] if candidates else source_dir

    args.output_dir.mkdir(parents=True, exist_ok=True)
    retained = []
    if source_dir.exists():
        for source in sorted(source_dir.glob("*.txt")):
            if source.name.startswith(SAFE_PREFIXES):
                shutil.copyfile(source, args.output_dir / source.name)
                retained.append(source.name)

    columns, rows = read_table(args.output_dir / "multiqc_general_stats.txt")
    links = "".join(
        f'<li><a href="{html.escape(args.output_dir.name)}/{html.escape(name)}">{html.escape(name)}</a></li>'
        for name in retained
    )
    document = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{html.escape(args.title)}</title>
<style>body{{font-family:system-ui,sans-serif;margin:2rem;color:#17202a}}table{{border-collapse:collapse;width:100%;font-size:.85rem}}th,td{{border:1px solid #ccd6dd;padding:.35rem;text-align:left}}th{{background:#e8f6f3}}li{{margin:.2rem 0}}.note{{color:#566573}}</style>
</head><body><h1>{html.escape(args.title)}</h1>
<p class="note">Deterministic offline summary of metrics parsed by MultiQC.</p>
<h2>General statistics</h2>{table_markup(columns, rows)}
<h2>Parsed data tables</h2><ul>{links}</ul>
</body></html>
"""
    args.report.write_text(document, encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
