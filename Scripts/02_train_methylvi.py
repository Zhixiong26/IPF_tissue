#!/usr/bin/env python3
"""Train MethylVI from the H5MU produced by 01_build_cov_methylvi_input.py."""
from __future__ import annotations
import argparse, json, os
from pathlib import Path
import anndata as ad
import mudata, numpy as np, scanpy as sc, scvi, torch
from scvi.external import METHYLVI

def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, default=Path(os.environ["IPF_MVI_INPUT"]))
    p.add_argument("--output", type=Path, default=Path(os.environ["IPF_MVI_RESULTS"]))
    p.add_argument("--batch-key", default=os.environ["IPF_BATCH_KEY"])
    p.add_argument("--epochs", type=int, default=500); p.add_argument("--batch-size", type=int, default=32)
    p.add_argument("--threads", type=int, default=int(os.environ["IPF_THREADS"])); p.add_argument("--seed", type=int, default=int(os.environ["IPF_SEED"]))
    p.add_argument("--accelerator", choices=("auto", "cpu", "gpu"), default="auto")
    args = p.parse_args(); args.output.mkdir(parents=True, exist_ok=True)
    scvi.settings.seed = args.seed; torch.set_num_threads(args.threads)
    mdata = mudata.read_h5mu(args.input); adata = mdata["mCG"]
    if args.batch_key not in adata.obs: raise KeyError(f"Batch key absent: {args.batch_key}")
    mc, cov = adata.layers["mc"], adata.layers["cov"]
    if np.any(mc > cov) or not np.issubdtype(mc.dtype, np.integer): raise ValueError("Invalid integer mc/cov layers")
    METHYLVI.setup_mudata(mdata, mc_layer="mc", cov_layer="cov", batch_key=args.batch_key,
                          methylation_contexts=["mCG"], modalities={"batch_key": "mCG"})
    model = METHYLVI(mdata, n_latent=20, n_hidden=128, n_layers=1, likelihood="betabinomial", dispersion="region")
    accelerator = "gpu" if args.accelerator == "auto" and torch.cuda.is_available() else args.accelerator
    if accelerator == "auto": accelerator = "cpu"
    model.train(max_epochs=args.epochs, early_stopping=True, batch_size=args.batch_size, accelerator=accelerator, devices=1)
    model.save(args.output / "model", overwrite=True, save_anndata=False)
    latent = np.asarray(model.get_latent_representation(batch_size=args.batch_size), dtype=np.float32)
    embedding = ad.AnnData(X=latent, obs=adata.obs.copy()); embedding.obsm["X_methylVI"] = latent
    sc.pp.neighbors(embedding, use_rep="X_methylVI", random_state=args.seed); sc.tl.umap(embedding, random_state=args.seed); sc.tl.leiden(embedding, key_added="methylVI_leiden", random_state=args.seed)
    embedding.write_h5ad(args.output / "methylvi_embedding.h5ad", compression="gzip")
    colors = [key for key in (args.batch_key, os.environ["IPF_CELLTYPE_KEY"], "methylVI_leiden") if key in embedding.obs]
    sc.pl.umap(embedding, color=colors, show=False, save=False, return_fig=True).savefig(args.output / "methylvi_umap.png", dpi=200, bbox_inches="tight")
    summary = {"input": str(args.input), "cells": embedding.n_obs, "features": adata.n_vars, "batch_key": args.batch_key, "accelerator": accelerator}
    (args.output / "run_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

if __name__ == "__main__": main()
