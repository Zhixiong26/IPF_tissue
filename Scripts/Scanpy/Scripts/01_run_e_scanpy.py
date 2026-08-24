#!/usr/bin/env python3
"""Run a reproducible Scanpy workflow on CYL/ZCP transcriptome E 10x matrices."""

import argparse
import json
import tempfile
import zipfile
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cyl-zip", type=Path, required=True)
    parser.add_argument("--zcp-zip", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--min-genes", type=int, default=200)
    parser.add_argument("--min-counts", type=int, default=500)
    parser.add_argument("--max-mt-percent", type=float, default=25.0)
    parser.add_argument("--min-cells-per-gene", type=int, default=3)
    parser.add_argument("--target-sum", type=float, default=10000.0)
    parser.add_argument("--n-hvg", type=int, default=3000)
    parser.add_argument("--n-pcs", type=int, default=30)
    parser.add_argument("--n-neighbors", type=int, default=15)
    parser.add_argument("--leiden-resolution", type=float, default=0.5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def read_10x_zip(path, cohort, temporary_root):
    if not path.is_file():
        raise FileNotFoundError(path)
    extract_root = Path(temporary_root) / cohort
    with zipfile.ZipFile(str(path)) as archive:
        archive.extractall(str(extract_root))
    matrix_dir = extract_root / "filtered_feature_bc_matrix"
    if not matrix_dir.is_dir():
        raise ValueError("10x matrix directory missing from %s" % path)
    adata = sc.read_10x_mtx(matrix_dir, var_names="gene_symbols", make_unique=True)
    adata.obs_names = pd.Index(["%s_%s" % (cohort, barcode) for barcode in adata.obs_names])
    adata.obs["cohort"] = cohort
    return adata


def save_current_figure(path):
    plt.gcf().savefig(str(path), dpi=180, bbox_inches="tight")
    plt.close("all")


def main():
    args = arguments()
    if (args.min_genes < 0 or args.min_counts < 0 or args.max_mt_percent < 0 or
            args.min_cells_per_gene < 1 or args.target_sum <= 0 or args.n_hvg < 1 or
            args.n_pcs < 1 or args.n_neighbors < 2):
        raise ValueError("QC thresholds and dimensionality parameters are invalid")

    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    h5ad_path = output_dir / "rna_e_scanpy.h5ad"
    if h5ad_path.exists() and not args.overwrite:
        raise FileExistsError("Refusing to overwrite %s; pass --overwrite" % h5ad_path)
    figures_dir = output_dir / "figures"
    figures_dir.mkdir(exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="scanpy_10x_", dir=str(output_dir)) as temporary_root:
        cyl = read_10x_zip(args.cyl_zip.resolve(), "CYL", temporary_root)
        zcp = read_10x_zip(args.zcp_zip.resolve(), "ZCP", temporary_root)
    input_cells = {"CYL": int(cyl.n_obs), "ZCP": int(zcp.n_obs)}
    adata = ad.concat([cyl, zcp], join="outer", merge="same", index_unique=None)
    if adata.obs_names.duplicated().any():
        raise ValueError("Duplicated sample-prefixed cell identifiers")
    adata.layers["counts"] = adata.X.copy()
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    adata.obs["pass_basic_qc"] = (
        (adata.obs["n_genes_by_counts"] >= args.min_genes) &
        (adata.obs["total_counts"] >= args.min_counts) &
        (adata.obs["pct_counts_mt"] <= args.max_mt_percent)
    )
    qc_columns = ["cohort", "n_genes_by_counts", "total_counts", "pct_counts_mt", "pass_basic_qc"]
    adata.obs[qc_columns].to_csv(output_dir / "cell_qc.tsv", sep="\t")
    sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
                 groupby="cohort", jitter=0.25, multi_panel=True, show=False)
    save_current_figure(figures_dir / "qc_violin_prefilter.png")
    sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", color="pct_counts_mt", show=False)
    save_current_figure(figures_dir / "qc_scatter_prefilter.png")

    adata = adata[adata.obs["pass_basic_qc"].to_numpy()].copy()
    if adata.n_obs < 3:
        raise RuntimeError("Fewer than three cells pass basic QC")
    sc.pp.filter_genes(adata, min_cells=args.min_cells_per_gene)
    if adata.n_vars < 2:
        raise RuntimeError("Fewer than two genes remain after filtering")
    sc.pp.normalize_total(adata, target_sum=args.target_sum)
    sc.pp.log1p(adata)
    adata.raw = adata
    hvg_count = min(args.n_hvg, adata.n_vars)
    sc.pp.highly_variable_genes(adata, n_top_genes=hvg_count, batch_key="cohort", flavor="seurat")
    if int(adata.var["highly_variable"].sum()) < 2:
        sc.pp.highly_variable_genes(adata, n_top_genes=hvg_count, flavor="seurat")
    adata = adata[:, adata.var["highly_variable"].to_numpy()].copy()

    n_comps = min(args.n_pcs, adata.n_obs - 1, adata.n_vars - 1)
    if n_comps < 2:
        raise RuntimeError("Too few cells or HVGs for PCA")
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=n_comps, svd_solver="arpack", random_state=args.seed)
    n_neighbors = min(args.n_neighbors, adata.n_obs - 1)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_comps, random_state=args.seed)
    sc.tl.leiden(adata, resolution=args.leiden_resolution, random_state=args.seed, key_added="leiden")
    sc.tl.umap(adata, random_state=args.seed)
    sc.pl.umap(adata, color=["cohort", "leiden"], show=False)
    save_current_figure(figures_dir / "umap_cohort_leiden.png")

    adata.uns["scanpy_run"] = {
        "inputs": {"CYL": str(args.cyl_zip.resolve()), "ZCP": str(args.zcp_zip.resolve())},
        "input_cells": input_cells,
        "min_genes": args.min_genes,
        "min_counts": args.min_counts,
        "max_mt_percent": args.max_mt_percent,
        "min_cells_per_gene": args.min_cells_per_gene,
        "target_sum": args.target_sum,
        "n_hvg_requested": args.n_hvg,
        "n_pcs": n_comps,
        "n_neighbors": n_neighbors,
        "leiden_resolution": args.leiden_resolution,
        "seed": args.seed,
        "scanpy_version": sc.__version__,
        "anndata_version": ad.__version__,
    }
    adata.write_h5ad(h5ad_path, compression="gzip")
    summary = {
        "input_cells": input_cells,
        "input_features": int(cyl.n_vars),
        "cells_after_qc": int(adata.n_obs),
        "hvg_after_selection": int(adata.n_vars),
        "cohort_cells_after_qc": adata.obs["cohort"].value_counts().sort_index().to_dict(),
        "leiden_clusters": int(adata.obs["leiden"].nunique()),
        "h5ad": str(h5ad_path),
    }
    with (output_dir / "run_summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print("Scanpy complete: %s" % h5ad_path, flush=True)
    print(json.dumps(summary, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
