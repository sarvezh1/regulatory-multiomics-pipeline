#!/usr/bin/env python3
import argparse
import json


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--route", choices=("raw", "processed"), required=True)
    parser.add_argument("--peaks", nargs="*", default=[])
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    peak_count = 0
    for path in args.peaks:
        with open(path, encoding="utf-8") as handle:
            peak_count += sum(1 for line in handle if line.strip() and not line.startswith("#"))

    if args.route == "raw":
        components = [
            {"component": "atac_processing", "status": "COMPLETE", "reason": "paired-end raw ATAC processing completed"},
            {
                "component": "atac_peak_calling",
                "status": "COMPLETE" if peak_count else "COMPLETE_EMPTY_RESULT",
                "reason": "MACS3 sample peak calling completed" if peak_count else "MACS3 completed with zero peaks",
            },
            {"component": "atac_count_contract", "status": "COMPLETE", "reason": "consensus regions and integer fragment matrix completed"},
        ]
    else:
        components = [
            {"component": "atac_processing", "status": "NOT_AVAILABLE_INPUT", "reason": "processed count-matrix route"},
            {"component": "atac_peak_calling", "status": "NOT_AVAILABLE_INPUT", "reason": "processed count-matrix route"},
            {"component": "atac_count_contract", "status": "COMPLETE", "reason": "external integer region matrix validated"},
        ]
    with open(args.output, "w", encoding="utf-8") as handle:
        json.dump({"components": components, "sample_peak_count": peak_count}, handle, sort_keys=True, indent=2)
        handle.write("\n")


if __name__ == "__main__":
    main()
