#!/usr/bin/env python3
"""Small standard-library helpers shared by the per-cell QC scripts."""

import csv
import gzip
import json
import os
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
