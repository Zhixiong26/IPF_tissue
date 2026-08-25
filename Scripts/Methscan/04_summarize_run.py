#!/usr/bin/env python3
"""Validate expected MethSCAn products and write a machine-readable run summary."""

import argparse
import csv
import gzip
import json
from pathlib import Path


def count_rows(path, delimiter=","):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(str(path), "rt", newline="") as handle:
        return max(sum(1 for _ in csv.reader(handle, delimiter=delimiter)) - 1, 0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    args = parser.parse_args()
    root = args.run_dir.resolve()
    required = {
        "manifest": root / "00_manifest/input_manifest.tsv",
        "scanpy_selected_manifest": root / "00_scanpy_selected/input_manifest.tsv",
        "scanpy_selection_summary": root / "00_scanpy_selected/scanpy_selection_summary.json",
        "prepared_stats": root / "01_prepared/cell_stats.csv",
        "filtered_stats": root / "02_filtered/cell_stats.csv",
        "vmrs": root / "03_scan/VMRs.bed",
        "methylated_sites": root / "04_matrix/methylated_sites.csv.gz",
        "total_sites": root / "04_matrix/total_sites.csv.gz",
        "methylation_fractions": root / "04_matrix/methylation_fractions.csv.gz",
        "mean_shrunken_residuals": root / "04_matrix/mean_shrunken_residuals.csv.gz",
        "embedding": root / "05_scanpy/tables/cell_embedding.tsv",
        "cell_missingness_qc": root / "05_scanpy/tables/cell_missingness_qc.tsv",
        "vmr_missingness_qc": root / "05_scanpy/tables/vmr_missingness_qc.tsv",
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
    with required["scanpy_selected_manifest"].open(newline="") as handle:
        selected_rows = list(csv.DictReader(handle, delimiter="\t"))
    if any(row.get("rna_cell_type") in ("", "NA") for row in selected_rows):
        raise ValueError("Scanpy-selected manifest contains missing or NA RNA cell types")
    rna_umap = root / "05_scanpy/figures/umap_rna_cell_type.png"
    if not rna_umap.is_file() or rna_umap.stat().st_size == 0:
        raise FileNotFoundError("RNA cell-type UMAP is missing/empty")
    required["rna_cell_type_umap"] = rna_umap
    input_cells = len(manifest_rows)
    scanpy_selected_cells = len(selected_rows)
    prepared_cells = count_rows(required["prepared_stats"])
    filtered_cells = count_rows(required["filtered_stats"])
    matrix_cells = count_rows(required["mean_shrunken_residuals"])
    scanpy_input_cells = int(scanpy["input_cells"])
    scanpy_cells = int(scanpy["retained_cells"])
    consistency_errors = []
    if prepared_cells != scanpy_selected_cells:
        consistency_errors.append("prepared_cells != scanpy_selected_cells")
    if filtered_cells > prepared_cells:
        consistency_errors.append("filtered_cells > prepared_cells")
    if matrix_cells != filtered_cells:
        consistency_errors.append("matrix_cells != filtered_cells")
    if scanpy_input_cells != matrix_cells:
        consistency_errors.append("scanpy input_cells != matrix_cells")
    if scanpy_cells > scanpy_input_cells:
        consistency_errors.append("scanpy retained_cells > scanpy input_cells")
    if consistency_errors:
        raise ValueError("Inconsistent MethSCAn stage counts: %s" % "; ".join(consistency_errors))
    summary = {
        "status": "complete", "run_dir": str(root), "input_cells": input_cells,
        "input_cells_by_sample": {}, "scanpy_selected_cells": scanpy_selected_cells,
        "prepared_cells": prepared_cells,
        "filtered_cells": filtered_cells, "matrix_cells": matrix_cells, "vmrs": n_vmrs,
        "scanpy_input_cells": scanpy_input_cells, "scanpy_cells": scanpy_cells,
        "scanpy_vmrs": scanpy["retained_regions"],
        "rna_annotation_available_cells": scanpy.get("rna_annotation_available_cells"),
        "filter_retention_fraction": filtered_cells / input_cells,
        "scanpy_retention_fraction": scanpy_cells / scanpy_input_cells,
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
        for key in ("status", "input_cells", "scanpy_selected_cells", "prepared_cells", "filtered_cells", "matrix_cells",
                    "vmrs", "scanpy_input_cells", "scanpy_cells", "scanpy_vmrs",
                    "rna_annotation_available_cells",
                    "filter_retention_fraction", "scanpy_retention_fraction"):
            handle.write("%s\t%s\n" % (key, summary[key]))
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
