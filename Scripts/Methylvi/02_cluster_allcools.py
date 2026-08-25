#!/usr/bin/env python3
"""Apply the yuanpei ALLCools mCG 5-kb clustering logic."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import warnings
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from ALLCools.clustering import (
    ConsensusClustering,
    binarize_matrix,
    filter_regions,
    lsi,
    significant_pc_test,
    tsne,
)
from ALLCools.mcds import MCDS


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mcds", type=Path, default=Path(os.environ["IPF_MCDS"]))
    parser.add_argument("--output", type=Path, default=Path(os.environ["IPF_ALLCOOLS_H5AD"]))
    parser.add_argument("--annotation", type=Path, default=Path(os.environ["IPF_ANNOTATION"]))
    parser.add_argument("--blacklist", type=Path, default=Path(os.environ["IPF_BLACKLIST"]))
    parser.add_argument("--blacklist-md5", default=os.environ["IPF_BLACKLIST_MD5"])
    parser.add_argument(
        "--blacklist-fraction", type=float,
        default=float(os.environ["IPF_BLACKLIST_FRACTION"]),
    )
    parser.add_argument(
        "--hypo-percent", type=float, default=float(os.environ["IPF_HYPO_PERCENT"]),
        help="Minimum percentage of cells with a binarized hypo-score for a bin to be retained",
    )
    parser.add_argument("--threads", type=int, default=int(os.environ["IPF_THREADS"]))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if not 0 < args.blacklist_fraction <= 1:
        raise ValueError("blacklist fraction must be in (0, 1]")
    if not 0 <= args.hypo_percent <= 100:
        raise ValueError("hypo percent must be in [0, 100]")
    if not args.blacklist.is_file():
        raise FileNotFoundError(args.blacklist)
    blacklist_md5 = hashlib.md5(args.blacklist.read_bytes()).hexdigest()
    if args.blacklist_md5 and blacklist_md5.lower() != args.blacklist_md5.lower():
        raise ValueError(
            f"Blacklist MD5 mismatch: expected={args.blacklist_md5} observed={blacklist_md5}"
        )
    warnings.filterwarnings("ignore", category=FutureWarning)
    mcds = MCDS.open(str(args.mcds), var_dim="chrom5k")
    bins_before_blacklist = int(mcds.get_index("chrom5k").size)
    mcds = mcds.remove_black_list_region(
        black_list_path=str(args.blacklist), f=args.blacklist_fraction,
    )
    bins_after_blacklist = int(mcds.get_index("chrom5k").size)
    print(
        f"Blacklist removed {bins_before_blacklist - bins_after_blacklist:,} / "
        f"{bins_before_blacklist:,} 5-kb bins (overlap >= {args.blacklist_fraction:g})",
        flush=True,
    )
    adata = mcds.get_score_adata(mc_type=os.environ["IPF_MC_CONTEXT"], quant_type="hypo-score")
    initial_shape = [int(adata.n_obs), int(adata.n_vars)]
    print(f"Initial matrix: {adata.n_obs:,} cells x {adata.n_vars:,} bins", flush=True)

    binarize_matrix(adata, cutoff=float(os.environ["IPF_BINARIZE_CUTOFF"]))
    filter_regions(adata, hypo_percent=args.hypo_percent)
    filtered_shape = [int(adata.n_obs), int(adata.n_vars)]
    print(f"Filtered matrix: {adata.n_obs:,} cells x {adata.n_vars:,} bins", flush=True)
    seed = int(os.environ["IPF_SEED"])
    requested_components = int(os.environ["IPF_LSI_COMPONENTS"])
    # scipy.sparse.linalg.svds (used by ARPACK) requires
    # 0 < k < min(matrix.shape). Keep the full-run value at 100 while making
    # small smoke datasets valid.
    lsi_components = min(requested_components, adata.n_obs - 1, adata.n_vars - 1)
    if lsi_components < 1:
        raise ValueError(
            f"LSI requires at least 2 cells and 2 retained bins; got {adata.shape}"
        )
    print(
        f"LSI components: {lsi_components} (requested {requested_components})",
        flush=True,
    )
    lsi(
        adata, n_components=lsi_components, algorithm="arpack",
        obsm="X_pca", random_state=seed,
    )
    n_components = significant_pc_test(
        adata, p_cutoff=float(os.environ["IPF_LSI_P_CUTOFF"]), update=True
    )
    print(f"Significant LSI components: {n_components}", flush=True)

    neighbors = int(os.environ["IPF_ALLCOOLS_NEIGHBORS"])
    sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=neighbors, random_state=seed)
    sc.tl.leiden(
        adata,
        resolution=float(os.environ["IPF_ALLCOOLS_LEIDEN_RESOLUTION"]),
        random_state=seed,
        key_added="leiden",
    )
    tsne(
        adata, obsm="X_pca", metric="euclidean", exaggeration=-1,
        perplexity=30, n_jobs=args.threads,
    )
    sc.tl.umap(adata, random_state=seed)
    consensus = ConsensusClustering(
        model=None, n_neighbors=neighbors, metric="euclidean", min_cluster_size=10,
        leiden_repeats=int(os.environ["IPF_CONSENSUS_LEIDEN_REPEATS"]),
        leiden_resolution=float(os.environ["IPF_CONSENSUS_LEIDEN_RESOLUTION"]),
        consensus_rate=0.5, random_state=seed, train_frac=0.5, train_max_n=500,
        max_iter=20, n_jobs=args.threads,
    )
    consensus.fit_predict(adata.obsm["X_pca"])
    adata.obs["L1"] = pd.Categorical(np.asarray(consensus.label).astype(str))
    adata.obs["L1_proba"] = np.asarray(consensus.label_proba, dtype=float)

    annotation = pd.read_csv(args.annotation, sep="\t", dtype=str)
    if annotation["cell_id"].duplicated().any():
        raise ValueError("Annotation contains duplicate cell_id values")
    annotation = annotation.set_index("cell_id")
    adata.obs["manual_celltype"] = (
        annotation.reindex(adata.obs_names)["manual_celltype"].fillna("Unknown").to_numpy()
    )
    adata.obs["cohort"] = adata.obs_names.to_series().str.split("_", n=1).str[0].to_numpy()
    adata.write_h5ad(args.output, compression="gzip")
    adata.obs.to_csv(args.output.parent / "cell_clusters.csv.gz")
    for basis in ("tsne", "umap"):
        sc.pl.embedding(
            adata, basis=basis,
            color=["L1", "L1_proba", "cohort", "manual_celltype"],
            show=False, wspace=0.35,
        )
        plt.savefig(args.output.parent / f"allcools_{basis}.png", dpi=300, bbox_inches="tight")
        plt.close("all")
    summary = {
        "initial_shape": initial_shape, "filtered_shape": filtered_shape,
        "blacklist": str(args.blacklist.resolve()), "blacklist_md5": blacklist_md5,
        "blacklist_fraction": args.blacklist_fraction,
        "bins_before_blacklist": bins_before_blacklist,
        "bins_after_blacklist": bins_after_blacklist,
        "blacklist_removed_bins": bins_before_blacklist - bins_after_blacklist,
        "hypo_percent": args.hypo_percent,
        "requested_lsi_components": requested_components,
        "computed_lsi_components": lsi_components,
        "significant_lsi_components": int(n_components), "neighbors": neighbors, "seed": seed,
    }
    (args.output.parent / "cluster_summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
