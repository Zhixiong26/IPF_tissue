#!/usr/bin/env python3
"""Small standard-library helpers shared by the per-cell QC scripts."""

import concurrent.futures as cf
import csv
import gzip
import json
import os
import time
from pathlib import Path


MISSING = {"", "NA", "NaN", "nan", "."}


def open_text(path):
    path = Path(path)
    return gzip.open(str(path), "rt") if path.name.endswith(".gz") else path.open("r")


def atomic_tsv(path, fieldnames, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(".%s.tmp.%s" % (path.name, os.getpid()))
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(".%s.tmp.%s" % (path.name, os.getpid()))
    with temporary.open("w") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
    temporary.replace(path)


def read_tsv(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def required_columns(rows, columns, path):
    if not rows:
        raise ValueError("Empty TSV: %s" % path)
    missing = set(columns) - set(rows[0])
    if missing:
        raise ValueError("%s is missing columns: %s" % (path, ", ".join(sorted(missing))))


def value_or_na(value):
    return "NA" if value is None else value


def parse_float(value):
    if value is None or str(value) in MISSING:
        return None
    return float(value)


def parse_int(value):
    if value is None or str(value) in MISSING:
        return None
    return int(value)


def cell_key(row):
    return row["sample_id"], row["barcode"]


def parallel_map_with_progress(function, tasks, workers, label, progress_interval):
    """Run tasks in input order while periodically reporting completed futures."""
    total = len(tasks)
    if not total:
        return []
    started = time.monotonic()
    print("%s progress: 0/%s completed; elapsed 0.0 min" % (label, total), flush=True)
    if workers == 1:
        results = []
        last_report = started
        for task in tasks:
            results.append(function(task))
            now = time.monotonic()
            if now - last_report >= progress_interval or len(results) == total:
                print("%s progress: %s/%s completed; elapsed %.1f min" %
                      (label, len(results), total, (now - started) / 60.0), flush=True)
                last_report = now
        return results

    results = [None] * total
    executor = cf.ProcessPoolExecutor(max_workers=workers)
    future_indexes = {
        executor.submit(function, task): index for index, task in enumerate(tasks)
    }
    pending = set(future_indexes)
    completed = 0
    next_report = started + progress_interval
    try:
        while pending:
            timeout = max(0.0, next_report - time.monotonic())
            done, pending = cf.wait(
                pending, timeout=timeout, return_when=cf.FIRST_COMPLETED
            )
            for future in done:
                results[future_indexes[future]] = future.result()
                completed += 1
            now = time.monotonic()
            if now >= next_report or not pending:
                print("%s progress: %s/%s completed; elapsed %.1f min" %
                      (label, completed, total, (now - started) / 60.0), flush=True)
                next_report = now + progress_interval
    finally:
        executor.shutdown()
    return results
