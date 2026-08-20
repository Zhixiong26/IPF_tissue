#!/usr/bin/env python3
"""Build a MethylVI H5MU directly from Bismark six-column .cov files.

The input files are *not* gzip compressed despite their ``allc.gz.cov``
suffix. Columns 5 and 6 are respectively methylated and unmethylated reads.
The script saves one sparse, resumable checkpoint per cell, selects informative
5-kb bins by across-cell methylation variance, then writes integer mc/cov
layers required by ``scvi.external.METHYLVI``.
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import gzip
import json
import os
import re
from collections import defaultdict
from pathlib import Path

import anndata as ad
import mudata
import numpy as np
import pandas as pd


def cell_id_from_path(path: Path) -> str:
    suffix = "_allc.gz.cov"
    if not path.name.endswith(suffix):
        raise ValueError(f"Unexpected coverage filename: {path.name}")
    return path.name[: -len(suffix)]


def bin_key(chrom: str, start_1based: int, bin_size: int) -> str:
    # Bismark .cov is 1-based inclusive. Store features as 0-based half-open.
    start = ((start_1based - 1) // bin_size) * bin_size
    return f"{chrom}:{start}-{start + bin_size}"


def build_checkpoint(task: tuple[str, str, str, int]) -> dict[str, object]:
    path_s, cell_id, checkpoint_s, bin_size = task
    path, checkpoint = Path(path_s), Path(checkpoint_s)
    if checkpoint.is_file() and checkpoint.stat().st_size:
        return {"cell_id": cell_id, "state": "reused"}
    sums: dict[str, list[int]] = {}
    with path.open("rt") as handle:
        for line_no, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 6:
                raise ValueError(f"{path}:{line_no}: expected 6 columns, got {len(fields)}")
            chrom, start_s, _end_s, _pct, mc_s, uc_s = fields
            mc, uc = int(mc_s), int(uc_s)
            if mc < 0 or uc < 0:
                raise ValueError(f"{path}:{line_no}: negative count")
            key = bin_key(chrom, int(start_s), bin_size)
            previous = sums.get(key)
            if previous is None:
                sums[key] = [mc, mc + uc]
            else:
                previous[0] += mc
                previous[1] += mc + uc
    features = np.asarray(sorted(sums), dtype=str)
    mc = np.asarray([sums[key][0] for key in features], dtype=np.uint32)
    cov = np.asarray([sums[key][1] for key in features], dtype=np.uint32)
    checkpoint.parent.mkdir(parents=True, exist_ok=True)
    temporary = checkpoint.with_suffix(".tmp.npz")
    np.savez_compressed(temporary, cell_id=cell_id, features=features, mc=mc, cov=cov)
    temporary.replace(checkpoint)
    return {"cell_id": cell_id, "state": "built", "features": len(features)}


def load_checkpoint(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with np.load(path, allow_pickle=False) as item:
        features, mc, cov = item["features"], item["mc"], item["cov"]
    if not (features.ndim == mc.ndim == cov.ndim == 1 and len(features) == len(mc) == len(cov)):
        raise ValueError(f"Invalid checkpoint: {path}")
    if np.any(mc > cov):
        raise ValueError(f"mc > cov in checkpoint: {path}")
    return features.astype(str), mc, cov


def feature_table(features: list[str]) -> pd.DataFrame:
    parsed = [re.fullmatch(r"(.+):(\d+)-(\d+)", value) for value in features]
    if any(item is None for item in parsed):
        raise ValueError("Malformed feature key")
    return pd.DataFrame(
        {"chrom": [x.group(1) for x in parsed], "start": [int(x.group(2)) for x in parsed],
         "end": [int(x.group(3)) for x in parsed]}, index=pd.Index(features, name="feature"),
    )


def read_chrom_sizes(path: Path) -> dict[str, int]:
    sizes: dict[str, int] = {}
    with path.open() as handle:
        for line_no, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 2 or int(fields[1]) <= 0:
                raise ValueError(f"{path}:{line_no}: expected chrom and positive size")
            sizes[fields[0]] = int(fields[1])
    if not sizes:
        raise ValueError(f"No chromosome sizes in {path}")
    return sizes


def is_canonical_feature(feature: str, chrom_sizes: dict[str, int]) -> bool:
    match = re.fullmatch(r"(.+):(\d+)-(\d+)", feature)
    return match is not None and match.group(1) in chrom_sizes and int(match.group(3)) <= chrom_sizes[match.group(1)]


def blacklist_overlaps(features: list[str], blacklist: Path, max_fraction: float) -> set[str]:
    """Return 5-kb bins whose ENCODE blacklist overlap exceeds max_fraction."""
    opener = gzip.open if blacklist.suffix == ".gz" else open
    intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with opener(blacklist, "rt") as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"{blacklist}:{line_no}: expected BED >=3 columns")
            intervals[fields[0]].append((int(fields[1]), int(fields[2])))
    for chrom, values in intervals.items():
        values.sort()
        merged: list[tuple[int, int]] = []
        for left, right in values:
            if merged and left <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(right, merged[-1][1]))
            else:
                merged.append((left, right))
        intervals[chrom] = merged
    bins_by_chrom: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    for feature in features:
        match = re.fullmatch(r"(.+):(\d+)-(\d+)", feature)
        if match is None:
            raise ValueError(f"Malformed feature key: {feature}")
        bins_by_chrom[match.group(1)].append((int(match.group(2)), int(match.group(3)), feature))
    rejected: set[str] = set()
    for chrom, bins in bins_by_chrom.items():
        bins.sort()
        chrom_intervals = intervals.get(chrom, [])
        interval_index = 0
        for start, end, feature in bins:
            while interval_index < len(chrom_intervals) and chrom_intervals[interval_index][1] <= start:
                interval_index += 1
            overlap, index = 0, interval_index
            while index < len(chrom_intervals) and chrom_intervals[index][0] < end:
                left, right = chrom_intervals[index]
                overlap += max(0, min(end, right) - max(start, left))
                index += 1
            if overlap / (end - start) > max_fraction:
                rejected.add(feature)
    return rejected


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cov-dir", type=Path, default=Path(os.environ["IPF_COV_DIR"]))
    parser.add_argument("--annotation", type=Path, default=Path(os.environ["IPF_ANNOTATION"]))
    parser.add_argument("--output", type=Path, default=Path(os.environ["IPF_MVI_INPUT"]))
    parser.add_argument("--work-dir", type=Path, default=Path(os.environ["IPF_MVI_ROOT"]) / "cov_checkpoints")
    parser.add_argument("--bin-size", type=int, default=int(os.environ["IPF_BIN_SIZE"]))
    parser.add_argument("--min-cells", type=int, default=int(os.environ["IPF_MIN_CELLS_PER_BIN"]))
    parser.add_argument("--min-cov", type=int, default=int(os.environ["IPF_MIN_COV_PER_BIN"]))
    parser.add_argument("--top-bins", type=int, default=int(os.environ["IPF_TOP_BINS"]))
    parser.add_argument("--threads", type=int, default=int(os.environ["IPF_THREADS"]))
    parser.add_argument("--chrom-sizes", type=Path, default=Path(os.environ["IPF_CHROM_SIZES"]))
    parser.add_argument("--blacklist", type=Path, default=Path(os.environ["IPF_BLACKLIST"]) if os.environ.get("IPF_BLACKLIST") else None)
    parser.add_argument("--blacklist-fraction", type=float, default=float(os.environ["IPF_BLACKLIST_FRACTION"]))
    parser.add_argument("--max-cells", type=int, help="Small test run only")
    parser.add_argument("--include-unannotated", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.bin_size <= 0 or args.top_bins <= 0 or args.min_cells <= 0 or args.min_cov <= 0:
        raise ValueError("bin-size, top-bins, min-cells and min-cov must be positive")
    if not 0 <= args.blacklist_fraction <= 1:
        raise ValueError("blacklist-fraction must be in [0, 1]")
    if not args.chrom_sizes.is_file():
        raise FileNotFoundError(args.chrom_sizes)
    if args.blacklist is not None and not args.blacklist.is_file():
        raise FileNotFoundError(args.blacklist)
    files = sorted(args.cov_dir.glob("*_allc.gz.cov"))
    if not files:
        raise FileNotFoundError(f"No *_allc.gz.cov files in {args.cov_dir}")
    annotation = pd.read_csv(args.annotation, sep="\t", dtype=str)
    required = {"cell_id", "manual_celltype"}
    if not required.issubset(annotation):
        raise ValueError(f"Annotation needs columns {sorted(required)}")
    annotation = annotation.drop_duplicates("cell_id").set_index("cell_id")
    cells = [cell_id_from_path(path) for path in files]
    keep = [cell in annotation.index and pd.notna(annotation.at[cell, "manual_celltype"]) and annotation.at[cell, "manual_celltype"] != "NA" for cell in cells]
    if not args.include_unannotated:
        files = [path for path, selected in zip(files, keep) if selected]
        cells = [cell for cell, selected in zip(cells, keep) if selected]
    if args.max_cells:
        files, cells = files[: args.max_cells], cells[: args.max_cells]
    if not cells:
        raise RuntimeError("No cells remain after annotation filtering")
    print(f"selected {len(cells):,} cells from {len(files):,} .cov files", flush=True)

    checkpoints = [args.work_dir / f"{cell}.npz" for cell in cells]
    tasks = [(str(path), cell, str(checkpoint), args.bin_size) for path, cell, checkpoint in zip(files, cells, checkpoints)]
    with cf.ProcessPoolExecutor(max_workers=args.threads) as pool:
        for index, result in enumerate(pool.map(build_checkpoint, tasks), 1):
            if index % 50 == 0 or index == len(tasks):
                print(f"checkpointed {index:,}/{len(tasks):,} ({result['state']})", flush=True)

    # Per-bin online moments of cell-level methylation ratios.  Bins observed
    # in few cells or with very low total coverage are excluded before ranking.
    stats: dict[str, list[float]] = defaultdict(lambda: [0.0, 0.0, 0.0, 0.0])
    for index, checkpoint in enumerate(checkpoints, 1):
        features, mc, cov = load_checkpoint(checkpoint)
        valid = cov >= args.min_cov
        for feature, ratio, depth in zip(features[valid], mc[valid] / cov[valid], cov[valid]):
            value = stats[feature]
            value[0] += 1  # observed cells
            value[1] += float(ratio)
            value[2] += float(ratio * ratio)
            value[3] += float(depth)
        if index % 200 == 0 or index == len(checkpoints):
            print(f"summarized {index:,}/{len(checkpoints):,} cells", flush=True)
    ranked = []
    for feature, (n, total, total_sq, depth) in stats.items():
        if n >= args.min_cells:
            variance = max(0.0, total_sq / n - (total / n) ** 2)
            ranked.append((variance, n, depth, feature))
    ranked.sort(reverse=True)
    chrom_sizes = read_chrom_sizes(args.chrom_sizes)
    before_canonical = len(ranked)
    ranked = [item for item in ranked if is_canonical_feature(item[3], chrom_sizes)]
    canonical_rejected = before_canonical - len(ranked)
    blacklist_rejected: set[str] = set()
    if args.blacklist is not None:
        blacklist_rejected = blacklist_overlaps([item[3] for item in ranked], args.blacklist, args.blacklist_fraction)
        print(f"blacklist removed {len(blacklist_rejected):,} candidate bins (>{args.blacklist_fraction:g} overlap)", flush=True)
    selected = [feature for _var, _n, _depth, feature in ranked if feature not in blacklist_rejected][: args.top_bins]
    if not selected:
        raise RuntimeError("No bins passed filtering; reduce --min-cells or --min-cov")
    selected.sort()
    feature_index = {feature: i for i, feature in enumerate(selected)}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.save(args.output.parent / "selected_5kb_bins.npy", np.asarray(selected, dtype=str))

    n_cells, n_features = len(cells), len(selected)
    max_count = 0
    for checkpoint in checkpoints:
        _features, _mc, cov = load_checkpoint(checkpoint)
        max_count = max(max_count, int(cov.max(initial=0)))
    dtype = np.uint16 if max_count <= np.iinfo(np.uint16).max else np.uint32
    mc_matrix = np.zeros((n_cells, n_features), dtype=dtype)
    cov_matrix = np.zeros((n_cells, n_features), dtype=dtype)
    for row, checkpoint in enumerate(checkpoints):
        features, mc, cov = load_checkpoint(checkpoint)
        indices = np.fromiter((feature_index.get(x, -1) for x in features), dtype=np.int64, count=len(features))
        keep = indices >= 0
        mc_matrix[row, indices[keep]] = mc[keep]
        cov_matrix[row, indices[keep]] = cov[keep]
    obs = pd.DataFrame(index=pd.Index(cells, name="cell_id"))
    obs["cohort"] = obs.index.to_series().str.split("_", n=1).str[0].to_numpy()
    obs["manual_celltype"] = annotation.reindex(cells)["manual_celltype"].fillna("Unknown").to_numpy()
    adata = ad.AnnData(X=None, obs=obs, var=feature_table(selected))
    adata.layers["mc"] = mc_matrix
    adata.layers["cov"] = cov_matrix
    temporary = args.output.with_suffix(".tmp.h5mu")
    mudata.MuData({"mCG": adata}).write_h5mu(temporary, compression="gzip")
    temporary.replace(args.output)
    summary = {"output": str(args.output), "cells": n_cells, "features": n_features, "bin_size": args.bin_size,
               "min_cells": args.min_cells, "min_cov": args.min_cov, "max_count": max_count, "dtype": np.dtype(dtype).name,
               "chrom_sizes": str(args.chrom_sizes), "noncanonical_candidate_bins_removed": canonical_rejected,
               "blacklist": str(args.blacklist) if args.blacklist else None, "blacklist_fraction": args.blacklist_fraction,
               "blacklist_candidate_bins_removed": len(blacklist_rejected),
               "cohorts": obs["cohort"].value_counts().to_dict(), "celltypes": obs["manual_celltype"].value_counts().to_dict()}
    (args.output.parent / "build_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
