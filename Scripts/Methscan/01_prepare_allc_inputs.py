#!/usr/bin/env python3
"""Discover, validate, select, and canonicalize ALLCools ALLC inputs."""

import argparse
import concurrent.futures as cf
import csv
import gzip
import json
import os
import re
import shutil
import tarfile
from collections import Counter
from pathlib import Path


ALLC_RE = re.compile(
    r"^(?P<sample>[^_]+)_(?P<barcode>[ACGT]{17})(?:\.allc\.tsv|_allc|\.allc)\.gz$"
)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--allc-source", type=Path, required=True,
                        help="ALLC directory or .tar/.tar.gz/.tgz archive")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--sample", action="append", dest="samples", required=True)
    parser.add_argument("--validation-records", type=int, default=10000,
                        help="Records checked per ALLC; 0 checks the complete file")
    parser.add_argument("--max-cells", type=int, default=0,
                        help="Balanced smoke subset; 0 retains all selected cells")
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--require-index", action="store_true",
                        help="Require a readable .tbi next to every source ALLC")
    parser.add_argument("--verify-only", action="store_true")
    return parser.parse_args()


def safe_extract(archive, destination):
    destination.mkdir(parents=True, exist_ok=False)
    root = destination.resolve()
    with tarfile.open(str(archive), "r:*") as handle:
        for member in handle.getmembers():
            target = (destination / member.name).resolve()
            if root != target and root not in target.parents:
                raise ValueError("Unsafe archive path: %s" % member.name)
            if member.issym() or member.islnk() or member.isdev():
                raise ValueError("Archive links/devices are not accepted: %s" % member.name)
        handle.extractall(str(destination))
    return destination


def archive_paths_below(directory):
    """Return supported archives staged below a project-local input directory."""
    return sorted(
        path for path in directory.rglob("*")
        if path.is_file() and (
            path.name.endswith(".tar")
            or path.name.endswith(".tar.gz")
            or path.name.endswith(".tgz")
        )
    )


def validate_allc(path, record_limit):
    records = 0
    contexts = Counter()
    previous = None
    with gzip.open(str(path), "rt") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                raise ValueError("%s:%d expected at least 6 columns" % (path, line_number))
            chrom, pos_s, strand, context, mc_s, cov_s = fields[:6]
            try:
                pos, mc, cov = int(pos_s), int(mc_s), int(cov_s)
            except ValueError as exc:
                raise ValueError("%s:%d invalid integer" % (path, line_number)) from exc
            if pos < 1 or strand not in ("+", "-") or mc < 0 or cov < mc:
                raise ValueError("%s:%d invalid coordinate/strand/count" % (path, line_number))
            key = (chrom, pos)
            if previous == key:
                raise ValueError("%s:%d duplicate coordinate" % (path, line_number))
            previous = key
            contexts[context] += 1
            records += 1
            if record_limit and records >= record_limit:
                break
    if records == 0:
        raise ValueError("Empty ALLC: %s" % path)
    return records, contexts


def validate_task(task):
    return validate_allc(task[0], task[1])


def balanced_subset(rows, max_cells):
    if not max_cells or len(rows) <= max_cells:
        return rows
    groups = {}
    for row in rows:
        groups.setdefault(row["sample_id"], []).append(row)
    selected = []
    while len(selected) < max_cells and any(groups.values()):
        for sample in sorted(groups):
            if groups[sample] and len(selected) < max_cells:
                selected.append(groups[sample].pop(0))
    return selected


def main():
    args = parse_args()
    if args.validation_records < 0 or args.max_cells < 0 or args.workers < 1:
        raise ValueError("validation-records/max-cells must be non-negative and workers positive")
    source = args.allc_source.resolve()
    if not source.exists():
        raise FileNotFoundError(source)
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise FileExistsError("ALLC intake output is not empty: %s" % args.output_dir)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    scan_roots = [source]
    source_kind = "directory"
    if source.is_file():
        source_kind = "archive"
        if args.verify_only:
            with tarfile.open(str(source), "r:*") as handle:
                names = handle.getnames()
            args.output_dir.rmdir()
            print(json.dumps({"status": "archive-readable", "members": len(names)}, indent=2))
            return
        scan_roots = [safe_extract(source, args.output_dir / "extracted")]
    elif not source.is_dir():
        raise ValueError("ALLC source is neither a directory nor an archive")
    else:
        archives = archive_paths_below(source)
        if archives:
            source_kind = "directory_with_archives"
            extracted_root = args.output_dir / "extracted"
            scan_roots.extend(
                safe_extract(archive, extracted_root / ("%03d_%s" % (index, archive.parent.name)))
                for index, archive in enumerate(archives, start=1)
            )

    samples = set(args.samples)
    candidates = []
    for scan_root in scan_roots:
        for path in sorted(scan_root.rglob("*.gz")):
            match = ALLC_RE.match(path.name)
            if match and match.group("sample") in samples:
                candidates.append((path, match.group("sample"), match.group("barcode")))
    if not candidates:
        raise ValueError("No matching ALLC files below %s" % source)

    selected_candidates = []
    seen = set()
    for path, sample, barcode in candidates:
        key = (sample, barcode)
        if key in seen:
            raise ValueError("Duplicate ALLC cell: %s_%s" % key)
        seen.add(key)
        selected_candidates.append({
            "sample_id": sample, "barcode": barcode, "cell_id": "%s_%s" % key,
            "source_path": str(path.resolve()),
        })
    selected_candidates = balanced_subset(selected_candidates, args.max_cells)
    validation_tasks = [
        (Path(row["source_path"]), args.validation_records) for row in selected_candidates
    ]
    with cf.ProcessPoolExecutor(max_workers=args.workers) as executor:
        validation_results = list(executor.map(validate_task, validation_tasks))

    rows = []
    context_counts = Counter()
    for candidate, validation in zip(selected_candidates, validation_results):
        path = Path(candidate["source_path"])
        sample, barcode = candidate["sample_id"], candidate["barcode"]
        index = Path(str(path) + ".tbi")
        if args.require_index and (not index.is_file() or not os.access(str(index), os.R_OK)):
            raise FileNotFoundError("Missing/read-protected ALLC index: %s" % index)
        checked, contexts = validation
        context_counts.update(contexts)
        cell_id = candidate["cell_id"]
        rows.append({
            "sample_id": sample, "barcode": barcode, "cell_id": cell_id,
            "source_path": str(path.resolve()), "source_index": str(index.resolve()) if index.exists() else "NA",
            "checked_records": checked,
        })
    if not rows:
        raise ValueError("No ALLC cells remain after selection")

    if args.verify_only:
        shutil.rmtree(str(args.output_dir))
        print(json.dumps({"status": "pass", "selected_total": len(rows)}, indent=2))
        return
    link_dir = args.output_dir / "input_links"
    link_dir.mkdir()
    for row in rows:
        link_path = link_dir / (row["cell_id"] + ".allc.gz")
        os.symlink(row["source_path"], str(link_path))
        row["link_path"] = str(link_path.absolute())
    fields = ["sample_id", "barcode", "cell_id", "source_path", "source_index",
              "link_path", "checked_records"]
    with (args.output_dir / "input_manifest.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    summary = {
        "status": "pass", "source": str(source), "source_kind": source_kind,
        "available_total": len(candidates), "selected_total": len(rows),
        "selected_by_sample": dict(sorted(Counter(row["sample_id"] for row in rows).items())),
        "validation_records_per_cell": args.validation_records,
        "validation_workers": args.workers,
        "observed_context_prefixes": dict(sorted(context_counts.items())),
    }
    with (args.output_dir / "manifest_summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
