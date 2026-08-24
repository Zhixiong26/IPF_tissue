#!/usr/bin/env python3
"""Count primary mapped BAM records and read1 records by CB tag."""

import argparse
import concurrent.futures as cf
import subprocess
from collections import Counter
from pathlib import Path

from qc_common import atomic_json, atomic_tsv


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--bam", action="append", default=[], help="Repeat for each BAM partition")
    parser.add_argument("--bam-dir", action="append", type=Path, default=[],
                        help="Repeat for BAM directories; files are selected with --bam-pattern")
    parser.add_argument("--bam-pattern", default="*_bismark_bt2_pe.bam")
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--cb-tag", default="CB")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--metrics", type=Path)
    parser.add_argument("--workers", type=int, default=16,
                        help="Parallel BAM partitions; keep modest because scanning is I/O-bound")
    parser.add_argument("--max-records", type=int, default=0, help="Testing only; 0 scans all records")
    return parser.parse_args()


def scan_bam(task):
    bam, samtools, cb_tag, max_records = task
    mapped_reads, mapped_pairs = Counter(), Counter()
    total = missing_cb = read2_records = 0
    command = [samtools, "view", "-F", "2308", str(bam)]
    process = subprocess.Popen(
        command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True
    )
    assert process.stdout is not None
    for line in process.stdout:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 11:
            process.kill()
            raise ValueError("Malformed SAM record from %s" % bam)
        flag = int(fields[1])
        barcode = None
        prefix = cb_tag + ":Z:"
        for tag in fields[11:]:
            if tag.startswith(prefix):
                barcode = tag[len(prefix):].upper()
                break
        total += 1
        if barcode is None:
            missing_cb += 1
        else:
            mapped_reads[barcode] += 1
            if flag & 64:
                mapped_pairs[barcode] += 1
            if flag & 128:
                read2_records += 1
        if max_records and total >= max_records:
            process.terminate()
            break
    stderr = process.stderr.read() if process.stderr else ""
    return_code = process.wait()
    if return_code and not (max_records and total >= max_records and return_code in (-15, 143)):
        raise RuntimeError("samtools failed for %s: %s" % (bam, stderr.strip()))
    return {"bam": str(bam), "mapped_reads": mapped_reads, "mapped_pairs": mapped_pairs,
            "total": total, "missing_cb": missing_cb, "read2_records": read2_records}


def main():
    args = arguments()
    if args.workers < 1 or args.max_records < 0:
        raise ValueError("--workers must be >0 and --max-records must be >=0")
    bam_files = [Path(path) for path in args.bam]
    for directory in args.bam_dir:
        if not directory.is_dir():
            raise NotADirectoryError(directory)
        bam_files.extend(sorted(directory.glob(args.bam_pattern)))
    bam_files = sorted(set(path.resolve() for path in bam_files))
    if not bam_files:
        raise ValueError("Provide --bam and/or --bam-dir with matching BAM files")
    mapped_reads, mapped_pairs = Counter(), Counter()
    total = missing_cb = read2_records = 0
    effective_workers = min(args.workers, len(bam_files))
    base_tasks = [(bam, args.samtools, args.cb_tag, 0) for bam in bam_files]
    if args.max_records:
        effective_workers = 1
        results = []
        for task in base_tasks:
            remaining = args.max_records - sum(item["total"] for item in results)
            if remaining <= 0:
                break
            results.append(scan_bam(task[:-1] + (remaining,)))
    elif effective_workers == 1:
        results = [scan_bam(task) for task in base_tasks]
    else:
        with cf.ProcessPoolExecutor(max_workers=effective_workers) as executor:
            results = list(executor.map(scan_bam, base_tasks))
    for result in results:
        mapped_reads.update(result["mapped_reads"])
        mapped_pairs.update(result["mapped_pairs"])
        total += result["total"]
        missing_cb += result["missing_cb"]
        read2_records += result["read2_records"]
        print("BAM complete: %s (%s primary mapped records)" %
              (result["bam"], result["total"]), flush=True)
    if total == 0:
        raise RuntimeError("No primary mapped records were read")
    if not mapped_reads:
        raise RuntimeError("No primary mapped record contains %s:Z" % args.cb_tag)
    rows = []
    for barcode in sorted(mapped_reads):
        rows.append({"sample_id": args.sample_id, "barcode": barcode,
                     "cell_id": "%s_%s" % (args.sample_id, barcode),
                     "mapped_pairs": mapped_pairs[barcode], "mapped_reads": mapped_reads[barcode]})
    atomic_tsv(args.output, ["sample_id", "barcode", "cell_id", "mapped_pairs", "mapped_reads"], rows)
    metrics = args.metrics or args.output.with_suffix(".metrics.json")
    atomic_json(metrics, {"sample_id": args.sample_id,
                          "bam_files": [str(path) for path in bam_files],
                          "primary_mapped_records_scanned": total,
                          "records_missing_cb": missing_cb,
                          "records_with_cb": total - missing_cb,
                          "read2_records_with_cb": read2_records,
                          "cells": len(mapped_reads), "cb_tag": args.cb_tag,
                          "workers_requested": args.workers,
                          "workers_effective": effective_workers,
                          "bam_files_scanned": len(results)})


if __name__ == "__main__":
    main()
