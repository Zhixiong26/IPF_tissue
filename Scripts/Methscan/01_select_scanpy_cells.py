#!/usr/bin/env python3
"""Select non-NA Scanpy-annotated ALLC cells before MethSCAn prepare."""

import argparse
import csv
import json
import os
from collections import Counter
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--scanpy-annotation", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--exclude-cell-type", default="NA")
    parser.add_argument("--max-cells", type=int, default=0,
                        help="Balanced smoke subset after Scanpy selection; 0 keeps all")
    return parser.parse_args()


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


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
    if args.max_cells < 0:
        raise ValueError("max-cells must be non-negative")
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise FileExistsError("Scanpy selection output is not empty: %s" % args.output_dir)

    manifest_rows = read_tsv(args.manifest)
    if not manifest_rows or "cell_id" not in manifest_rows[0]:
        raise ValueError("Manifest must contain at least one cell_id row")
    if len({row["cell_id"] for row in manifest_rows}) != len(manifest_rows):
        raise ValueError("Manifest cell_id values must be unique")

    annotation_rows = read_tsv(args.scanpy_annotation)
    required = {"cell_id", "cell_type"}
    if not annotation_rows or not required.issubset(annotation_rows[0]):
        raise ValueError("Scanpy annotation must contain cell_id and cell_type")
    if len({row["cell_id"] for row in annotation_rows}) != len(annotation_rows):
        raise ValueError("Scanpy annotation cell_id values must be unique")
    annotation = {row["cell_id"]: row for row in annotation_rows}

    selected, excluded = [], []
    for row in manifest_rows:
        annotation_row = annotation.get(row["cell_id"])
        if annotation_row is None:
            excluded.append({"cell_id": row["cell_id"], "sample_id": row["sample_id"],
                             "reason": "not_in_scanpy"})
            continue
        if annotation_row["cell_type"] == args.exclude_cell_type:
            excluded.append({"cell_id": row["cell_id"], "sample_id": row["sample_id"],
                             "reason": "excluded_scanpy_cell_type_%s" % args.exclude_cell_type})
            continue
        selected.append(dict(row, **{
            "rna_cohort": annotation_row.get("cohort", ""),
            "rna_sample": annotation_row.get("sample", ""),
            "rna_leiden": annotation_row.get("leiden", ""),
            "rna_cell_type": annotation_row["cell_type"],
            "rna_annotation_available": "True",
        }))
    selected = balanced_subset(selected, args.max_cells)
    if not selected:
        raise ValueError("No ALLC cells remain after Scanpy non-NA selection")

    args.output_dir.mkdir(parents=True)
    links = args.output_dir / "input_links"
    links.mkdir()
    fields = list(manifest_rows[0].keys()) + [
        "rna_cohort", "rna_sample", "rna_leiden", "rna_cell_type", "rna_annotation_available"
    ]
    for row in selected:
        source = Path(row["source_path"])
        if not source.is_file():
            raise FileNotFoundError(source)
        link = links / (row["cell_id"] + ".allc.gz")
        os.symlink(str(source), str(link))
        row["link_path"] = str(link.absolute())
    with (args.output_dir / "input_manifest.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(selected)
    with (args.output_dir / "allc_excluded_by_scanpy.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["cell_id", "sample_id", "reason"], delimiter="\t")
        writer.writeheader()
        writer.writerows(excluded)
    summary = {
        "input_allc_cells": len(manifest_rows),
        "scanpy_selected_cells": len(selected),
        "selected_by_sample": dict(sorted(Counter(row["sample_id"] for row in selected).items())),
        "excluded_by_reason": dict(sorted(Counter(row["reason"] for row in excluded).items())),
        "excluded_scanpy_cell_type": args.exclude_cell_type,
    }
    with (args.output_dir / "scanpy_selection_summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
