#!/usr/bin/env python3
"""Validate expected MethSCAn products and write a machine-readable run summary."""

import argparse
import csv
import json
from pathlib import Path


def count_rows(path, delimiter=","):
    with path.open(newline="") as handle:
        return max(sum(1 for _ in csv.reader(handle, delimiter=delimiter)) - 1, 0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    args = parser.parse_args()
    root = args.run_dir.resolve()
    required = {
        "manifest": root / "00_manifest/input_manifest.tsv",
        "prepared_stats": root / "01_prepared/cell_stats.csv",
        "filtered_stats": root / "02_filtered/cell_stats.csv",
        "vmrs": root / "03_scan/VMRs.bed",
        "methylated_sites": root / "04_matrix/methylated_sites.csv.gz",
        "total_sites": root / "04_matrix/total_sites.csv.gz",
        "methylation_fractions": root / "04_matrix/methylation_fractions.csv.gz",
        "mean_shrunken_residuals": root / "04_matrix/mean_shrunken_residuals.csv.gz",
        "embedding": root / "05_scanpy/tables/cell_embedding.tsv",
        "h5ad": root / "05_scanpy/objects/methscan_vmr.h5ad",
        "scanpy_parameters": root / "05_scanpy/run_parameters.json",
    }
    missing = [name for name, path in required.items() if not path.is_file() or path.stat().st_size == 0]
    if missing:
        raise FileNotFoundError("Missing/empty outputs: %s" % ", ".join(missing))
    with required["manifest"].open(newline="") as handle:
        manifest_rows = list(csv.DictReader(handle, delimiter="\t"))
    with required["vmrs"].open() as handle:
        n_vmrs = sum(1 for line in handle if line.strip() and not line.startswith("#"))
    with required["scanpy_parameters"].open() as handle:
        scanpy = json.load(handle)
    summary = {
        "status": "complete", "run_dir": str(root), "input_cells": len(manifest_rows),
        "input_cells_by_sample": {}, "prepared_cells": count_rows(required["prepared_stats"]),
        "filtered_cells": count_rows(required["filtered_stats"]), "vmrs": n_vmrs,
        "scanpy_cells": scanpy["retained_cells"], "scanpy_vmrs": scanpy["retained_regions"],
        "outputs": {name: str(path) for name, path in required.items()},
    }
    for row in manifest_rows:
        sample = row["sample_id"]
        summary["input_cells_by_sample"][sample] = summary["input_cells_by_sample"].get(sample, 0) + 1
    with (root / "run_summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    with (root / "run_summary.tsv").open("w") as handle:
        handle.write("metric\tvalue\n")
        for key in ("status", "input_cells", "prepared_cells", "filtered_cells", "vmrs", "scanpy_cells", "scanpy_vmrs"):
            handle.write("%s\t%s\n" % (key, summary[key]))
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
