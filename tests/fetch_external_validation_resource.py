#!/usr/bin/env python3
"""Fetch a pinned validation resource with resumable HTTP byte ranges."""

import argparse
import concurrent.futures
import hashlib
import os
import urllib.request
from pathlib import Path


def download_part(url, part, start, end):
    part.parent.mkdir(parents=True, exist_ok=True)
    existing = part.stat().st_size if part.exists() else 0
    expected = end - start + 1
    if existing > expected:
        raise ValueError(f"oversized partial download: {part}")
    if existing == expected:
        return
    request = urllib.request.Request(
        url,
        headers={
            "Range": f"bytes={start + existing}-{end}",
        },
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        if response.status != 206:
            raise RuntimeError(f"server did not honor byte range for {part.name}")
        with part.open("ab") as handle:
            while True:
                chunk = response.read(1024 * 1024)
                if not chunk:
                    break
                handle.write(chunk)
    if part.stat().st_size != expected:
        raise RuntimeError(f"partial download has unexpected size: {part}")


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--url", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--size", required=True, type=int)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    if args.size <= 0 or args.workers <= 0:
        raise ValueError("size and workers must be positive")

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    chunk_size = (args.size + args.workers - 1) // args.workers
    ranges = []
    for index in range(args.workers):
        start = index * chunk_size
        if start >= args.size:
            break
        end = min(args.size - 1, start + chunk_size - 1)
        part = output.with_name(f"{output.name}.part.{index:03d}")
        ranges.append((args.url, part, start, end))

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = [executor.submit(download_part, *item) for item in ranges]
        for future in concurrent.futures.as_completed(futures):
            future.result()

    temporary = output.with_name(f"{output.name}.assembling")
    with temporary.open("wb") as destination:
        for _, part, _, _ in ranges:
            with part.open("rb") as source:
                while True:
                    block = source.read(1024 * 1024)
                    if not block:
                        break
                    destination.write(block)
    if temporary.stat().st_size != args.size:
        raise RuntimeError("assembled download has unexpected size")
    os.replace(temporary, output)
    for _, part, _, _ in ranges:
        part.unlink()
    print(f"path={output.as_posix()} bytes={output.stat().st_size} sha256={sha256(output)}")


if __name__ == "__main__":
    main()
