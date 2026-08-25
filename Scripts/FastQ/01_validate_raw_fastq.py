#!/usr/bin/env python3
"""Validate linked paired-end FASTQ layout without scanning entire gzip streams."""
import argparse
import csv
import gzip
import json
import os
import re
from collections import Counter, defaultdict
from pathlib import Path


FASTQ_RE = re.compile(
    r"^(?P<pair_id>.+)_R(?P<read>[12])_001\.(?:fastq|fq)\.gz$",
    re.IGNORECASE,
)


def first_record(path):
    with gzip.open(path, "rt", encoding="ascii", errors="strict") as handle:
        header = handle.readline().rstrip("\r\n")
        sequence = handle.readline().rstrip("\r\n")
        plus = handle.readline().rstrip("\r\n")
        quality = handle.readline().rstrip("\r\n")
    if not header.startswith("@"):
        raise ValueError("first FASTQ header does not start with '@'")
    if not sequence:
        raise ValueError("first FASTQ sequence is empty")
    if not plus.startswith("+"):
        raise ValueError("first FASTQ separator does not start with '+'")
    if len(sequence) != len(quality):
        raise ValueError("first sequence and quality lengths differ")
    read_name = header[1:].split(maxsplit=1)[0]
    read_name = re.sub(r"/[12]$", "", read_name)
    return read_name, len(sequence)


def find_fastqs(root):
    paths = []
    for directory, _subdirs, filenames in os.walk(root, followlinks=True):
        for filename in filenames:
            if FASTQ_RE.match(filename):
                paths.append(Path(directory) / filename)
    return sorted(paths)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--raw-dir",
        type=Path,
        default=Path("/home/lijia/luozhixiong/IPF_tissue/Data/Raw_fastq"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("/home/lijia/luozhixiong/IPF_tissue/Results/FastQ/raw_fastq_validation"),
    )
    args = parser.parse_args()

    raw_dir = args.raw_dir.resolve(strict=True)
    if not raw_dir.is_dir() or not os.access(raw_dir, os.R_OK | os.X_OK):
        raise PermissionError(f"Raw FASTQ directory is not readable: {raw_dir}")

    top_links = sorted(path for group in args.raw_dir.iterdir() if group.is_dir() for path in group.iterdir())
    broken_links = [str(path) for path in top_links if path.is_symlink() and not path.exists()]
    files = find_fastqs(args.raw_dir)
    if not files:
        raise FileNotFoundError(f"No paired-end *.fastq.gz or *.fq.gz files under {args.raw_dir}")

    rows = []
    pair_members = defaultdict(dict)
    errors = []
    total_bytes = 0

    for path in files:
        relative = path.relative_to(args.raw_dir)
        match = FASTQ_RE.match(path.name)
        assert match is not None
        read = match.group("read")
        pair_id = match.group("pair_id")
        size = path.stat().st_size
        total_bytes += size
        status = "ok"
        read_name = ""
        read_length = ""
        try:
            if size <= 0:
                raise ValueError("file is empty")
            if not os.access(path, os.R_OK):
                raise PermissionError("file is not readable")
            read_name, read_length = first_record(path)
        except Exception as exc:  # keep a complete manifest even when one file fails
            status = "error"
            errors.append(f"{relative}: {exc}")

        parts = relative.parts
        group = parts[0] if len(parts) >= 2 else ""
        batch = parts[1] if len(parts) >= 3 else ""
        subset = "/".join(parts[2:-1])
        row = {
            "group": group,
            "batch": batch,
            "subset": subset,
            "pair_id": pair_id,
            "read": read,
            "filename": path.name,
            "path": str(path),
            "bytes": size,
            "first_read_name": read_name,
            "first_read_length": read_length,
            "status": status,
        }
        rows.append(row)
        key = (str(path.parent), pair_id)
        if read in pair_members[key]:
            errors.append(f"duplicate R{read} for pair {pair_id} in {path.parent}")
        pair_members[key][read] = row

    for (parent, pair_id), members in sorted(pair_members.items()):
        if set(members) != {"1", "2"}:
            errors.append(f"incomplete pair {pair_id} in {parent}: reads={sorted(members)}")
            continue
        r1, r2 = members["1"], members["2"]
        if r1["first_read_name"] != r2["first_read_name"]:
            errors.append(f"first R1/R2 read names differ for {pair_id} in {parent}")

    if broken_links:
        errors.extend(f"broken top-level link: {path}" for path in broken_links)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = args.output_dir / "raw_fastq_manifest.tsv"
    with manifest_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    group_counts = Counter(str(row["group"]) for row in rows)
    group_bytes = Counter()
    for row in rows:
        group_bytes[str(row["group"])] += int(row["bytes"])
    summary = {
        "raw_dir": str(args.raw_dir),
        "resolved_raw_dir": str(raw_dir),
        "top_level_entries": len(top_links),
        "broken_top_level_links": broken_links,
        "fastq_files": len(rows),
        "candidate_pairs": len(pair_members),
        "files_by_group": dict(sorted(group_counts.items())),
        "bytes_by_group": dict(sorted(group_bytes.items())),
        "total_bytes": total_bytes,
        "validation_scope": "layout, size/readability, gzip decoding of first record, FASTQ syntax, and first R1/R2 read-name match; not full-stream gzip integrity",
        "status": "pass" if not errors else "fail",
        "errors": errors,
        "manifest": str(manifest_path),
    }
    summary_path = args.output_dir / "raw_fastq_validation_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
