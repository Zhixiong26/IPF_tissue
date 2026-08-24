#!/usr/bin/env python3
"""Count paired-end FASTQ fragments per cell from the R1 cell barcode."""

import argparse
import csv
import re
from collections import Counter
from pathlib import Path

from qc_common import atomic_json, atomic_tsv, open_text, parallel_map_with_progress


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--r1", action="append", required=True, help="Repeat for each lane/partition")
    parser.add_argument("--r2", action="append", required=True, help="Repeat in the same order as --r1")
    parser.add_argument("--barcode-source", choices=("header", "r1", "r2"), default="r1",
                        help="Project default: r1")
    parser.add_argument("--barcode-regex", help="Regex with exactly one capture group; required for header mode")
    parser.add_argument("--barcode-start", type=int, default=0,
                        help="Zero-based sequence start; project default: 0")
    parser.add_argument("--barcode-length", type=int, default=17,
                        help="Project default: 17")
    parser.add_argument(
        "--whitelist", type=Path,
        help="One barcode per line, or a step-02 TSV with sample_id/barcode columns",
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--metrics", type=Path)
    parser.add_argument("--workers", type=int, default=2,
                        help="Parallel FASTQ partitions; project maximum is normally 2 per sample")
    parser.add_argument("--progress-interval", type=int, default=180, metavar="SECONDS",
                        help="Heartbeat interval for parallel progress logs (default: 180)")
    parser.add_argument("--max-pairs", type=int, default=0, help="Testing only; 0 scans all pairs")
    return parser.parse_args()


def records(path):
    with open_text(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            sequence, plus, quality = handle.readline(), handle.readline(), handle.readline()
            if not sequence or not plus or not quality:
                raise ValueError("Truncated FASTQ record in %s" % path)
            if not header.startswith("@") or not plus.startswith("+"):
                raise ValueError("Invalid FASTQ record in %s: %s" % (path, header.rstrip()))
            if len(sequence.rstrip()) != len(quality.rstrip()):
                raise ValueError("Sequence/quality length mismatch in %s: %s" % (path, header.rstrip()))
            yield header.rstrip(), sequence.rstrip()


def read_name(header):
    return header[1:].split()[0]


def load_whitelist(path, sample_id, barcode_length):
    if path is None:
        return None
    with path.open(newline="") as handle:
        first = handle.readline()
        handle.seek(0)
        header = first.rstrip("\n").split("\t")
        if "barcode" in header:
            reader = csv.DictReader(handle, delimiter="\t")
            values = {
                row["barcode"].upper() for row in reader
                if row.get("barcode") and (not row.get("sample_id") or row["sample_id"] == sample_id)
            }
        else:
            values = {
                line.strip().split()[0].upper() for line in handle
                if line.strip() and not line.startswith("#")
            }
    if not values:
        raise ValueError("No whitelist barcode for sample %s in %s" % (sample_id, path))
    invalid = sorted(value for value in values if not re.fullmatch(r"[ACGT]{%d}" % barcode_length, value))
    if invalid:
        raise ValueError("Invalid whitelist barcode(s), first examples: %s" % invalid[:5])
    return values


def count_fastq_pair(task):
    (r1_path, r2_path, barcode_source, barcode_regex, barcode_start,
     barcode_length, whitelist, max_pairs) = task
    pattern = re.compile(barcode_regex) if barcode_source == "header" else None
    counts = Counter()
    total = invalid_barcode = not_in_whitelist = 0
    examples = []
    left, right = records(r1_path), records(r2_path)
    while True:
        try:
            one = next(left)
        except StopIteration:
            one = None
        try:
            two = next(right)
        except StopIteration:
            two = None
        if one is None or two is None:
            if one is not None or two is not None:
                raise ValueError("R1/R2 record counts differ: %s, %s" % (r1_path, r2_path))
            break
        if read_name(one[0]) != read_name(two[0]):
            raise ValueError("R1/R2 name mismatch: %s != %s" % (one[0], two[0]))
        if barcode_source == "header":
            match = pattern.search(one[0])
            barcode = match.group(1).upper() if match else None
        else:
            sequence = one[1] if barcode_source == "r1" else two[1]
            end = barcode_start + barcode_length
            barcode = sequence[barcode_start:end].upper() if len(sequence) >= end else None
        total += 1
        if barcode is None or not re.fullmatch(r"[ACGT]{%d}" % barcode_length, barcode):
            invalid_barcode += 1
            if len(examples) < 5:
                examples.append(one[0])
        elif whitelist and barcode not in whitelist:
            not_in_whitelist += 1
        else:
            counts[barcode] += 1
        if max_pairs and total >= max_pairs:
            break
    return {"r1": r1_path, "r2": r2_path, "counts": counts, "total": total,
            "invalid_barcode": invalid_barcode, "not_in_whitelist": not_in_whitelist,
            "examples": examples}


def main():
    args = arguments()
    if len(args.r1) != len(args.r2):
        raise ValueError("--r1 and --r2 counts differ")
    if (args.max_pairs < 0 or args.barcode_length < 1 or args.workers < 1 or
            args.progress_interval < 1):
        raise ValueError("--max-pairs must be >=0; --barcode-length/--workers/"
                         "--progress-interval must be >0")
    pattern = None
    if args.barcode_source == "header":
        if not args.barcode_regex:
            raise ValueError("--barcode-regex is required for header mode")
        pattern = re.compile(args.barcode_regex)
        if pattern.groups != 1:
            raise ValueError("--barcode-regex must contain exactly one capture group")
    elif args.barcode_start is None or args.barcode_start < 0:
        raise ValueError("--barcode-start >= 0 is required for r1/r2 sequence mode")

    whitelist = load_whitelist(args.whitelist, args.sample_id, args.barcode_length)
    counts = Counter()
    total = invalid_barcode = not_in_whitelist = 0
    examples = []
    base_tasks = [(r1, r2, args.barcode_source, args.barcode_regex, args.barcode_start,
                   args.barcode_length, whitelist, 0) for r1, r2 in zip(args.r1, args.r2)]
    effective_workers = min(args.workers, len(base_tasks))
    if args.max_pairs:
        effective_workers = 1
        results = []
        for task in base_tasks:
            remaining = args.max_pairs - sum(item["total"] for item in results)
            if remaining <= 0:
                break
            results.append(count_fastq_pair(task[:-1] + (remaining,)))
    else:
        results = parallel_map_with_progress(
            count_fastq_pair, base_tasks, effective_workers, "FASTQ", args.progress_interval
        )
    for result in results:
        counts.update(result["counts"])
        total += result["total"]
        invalid_barcode += result["invalid_barcode"]
        not_in_whitelist += result["not_in_whitelist"]
        examples.extend(result["examples"][:max(0, 5 - len(examples))])
        print("FASTQ complete: %s (%s pairs)" % (result["r1"], result["total"]), flush=True)
    if total == 0:
        raise RuntimeError("No FASTQ pairs were read")
    if not counts and not whitelist:
        raise RuntimeError("No cell barcode matched the explicit rule; examples: %s" % examples)
    if whitelist:
        for barcode in whitelist:
            counts.setdefault(barcode, 0)
    rows = [
        {"sample_id": args.sample_id, "barcode": barcode,
         "cell_id": "%s_%s" % (args.sample_id, barcode),
         "input_pairs": count, "input_reads": count * 2}
        for barcode, count in sorted(counts.items())
    ]
    atomic_tsv(args.output, ["sample_id", "barcode", "cell_id", "input_pairs", "input_reads"], rows)
    metrics = args.metrics or args.output.with_suffix(".metrics.json")
    assigned = sum(counts.values())
    atomic_json(metrics, {"sample_id": args.sample_id, "pairs_scanned": total,
                          "pairs_assigned": assigned,
                          "pairs_invalid_barcode": invalid_barcode,
                          "pairs_not_in_whitelist": not_in_whitelist,
                          "cells_in_output": len(counts),
                          "cells_with_input_pairs": sum(value > 0 for value in counts.values()),
                          "whitelist": str(args.whitelist.resolve()) if args.whitelist else None,
                          "barcode_source": args.barcode_source,
                          "barcode_regex": args.barcode_regex,
                          "barcode_start": args.barcode_start,
                          "barcode_length": args.barcode_length,
                          "workers_requested": args.workers,
                          "workers_effective": effective_workers,
                          "progress_interval_seconds": args.progress_interval,
                          "fastq_partitions": len(results)})


if __name__ == "__main__":
    main()
