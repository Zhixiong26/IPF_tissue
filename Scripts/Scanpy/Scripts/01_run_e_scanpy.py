#!/usr/bin/env python3
"""Run the confirmed CYL/ZCP transcriptome-E Scanpy workflow."""

from __future__ import annotations

import argparse
import json
import platform
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
import scanpy.external as sce
import scrublet as scr


PROJECT_DIR = Path("/home/lijia/luozhixiong/IPF_tissue")
DEFAULT_INPUTS = {
    "CYL": PROJECT_DIR / "Data/Matrix/25100718_CYL_E/filtered_feature_bc_matrix.zip",
    "ZCP": PROJECT_DIR / "Data/Matrix/25100718_ZCP_E/filtered_feature_bc_matrix.zip",
}

MARKER_GENES = {
    "Pan_epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19"],
    "AT2": ["SFTPC", "SFTPB", "SFTPA1", "SFTPA2", "ABCA3", "LPCAT1"],
    "AT1": ["AGER", "CAV1", "CAV2", "PDPN", "HOPX", "EMP2", "AQP5"],
    "Secretory": ["SCGB1A1", "SCGB3A1", "SCGB3A2", "BPIFB1", "MUC4", "WFDC2", "SLPI"],
    "Basal": ["KRT5", "KRT14", "KRT15", "KRT17", "TP63", "MIR205HG"],
    "Ciliated": ["FOXJ1", "PIFO", "TPPP3", "DNAH11", "CFAP46", "HYDIN"],
    "Macrophage": ["LST1", "TYROBP", "FCER1G", "CD68", "C1QA", "C1QB", "C1QC", "MRC1", "CD163", "PPARG"],
    "Fibroblast": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA", "COL6A3"],
    "Pan_endothelial": ["PECAM1", "VWF", "KDR", "EMCN", "ENG", "ESAM", "RAMP2", "PLVAP"],
    "Capillary_endothelial": ["CA4", "RGCC", "EMCN", "GPIHBP1", "BTNL9", "EDNRB", "EPAS1"],
    "Lymphatic_endothelial": ["PROX1", "PDPN", "LYVE1", "FLT4", "CCL21", "MMRN1", "RELN"],
    "Smooth_muscle_mural": ["ACTA2", "TAGLN", "MYH11", "LMOD1", "CNN1", "CARMN", "PDGFRB", "RGS5"],
    "T_cell": ["CD3D", "CD3E", "TRAC", "BCL11B", "ITK", "CD247", "IL7R"],
    "Mast": ["KIT", "CPA3", "TPSAB1", "TPSB2", "HDC", "MS4A2", "HPGDS"],
    "Cycling": ["MKI67", "TOP2A", "UBE2C", "CENPF", "BIRC5", "RRM2", "ANLN", "ECT2", "DIAPH3"],
}

DOTPLOT_MARKERS = {
    "AT2": ["SFTPC", "ABCA3", "NAPSA", "LPCAT1"],
    "Secretory epithelial": ["NEDD4L", "SFTA3", "SCNN1B", "GPRC5A"],
    "Fibroblasts": ["COL1A1", "COL1A2", "DCN", "COL3A1"],
    "Ciliated cells": ["FOXJ1", "DNAH11", "PIFO", "CFAP46"],
    "Secretory / mucous epithelial": ["BPIFB1", "MUC4", "WFDC2", "TMC5"],
    "Macrophages": ["LST1", "C1QA", "MRC1", "CD163"],
    "AT1-like": ["AGER", "CAV1", "HOPX", "AQP5"],
    "Endothelial cells": ["PECAM1", "VWF", "KDR", "EMCN"],
    "Basal cells": ["KRT5", "KRT15", "TP63", "KRT17"],
    "Smooth muscle / mural cells": ["ACTA2", "MYH11", "CARMN", "PDGFRB"],
    "MT-high AT2-like": ["MT-CO1", "MT-ND1", "MT-ND4", "MT-CYB"],
    "T cells": ["CD3D", "CD3E", "TRAC", "BCL11B"],
    "Cycling cells": ["MKI67", "TOP2A", "RRM2", "ANLN"],
    "Lymphatic endothelial cells": ["PROX1", "FLT4", "CCL21", "LYVE1"],
    "NA": ["EPCAM", "KRT8", "COL11A1", "CEMIP"],
}

CELL_TYPE_ORDER = list(DOTPLOT_MARKERS)

CLUSTER_ANNOTATIONS = {
    "0": ("AT2", "high", "SFTPC/SFTPB/ABCA3/LPCAT1"),
    "1": ("Secretory epithelial", "medium", "NEDD4L/SFTPB/SFTA3 with epithelial localization"),
    "2": ("Fibroblasts", "high", "COL1A2/COL3A1/COL5A1/COL6A3/PDGFRA"),
    "3": ("Ciliated cells", "high", "CFAP/DNAH/HYDIN"),
    "4": ("Secretory / mucous epithelial", "high", "BPIFB1/MUC4/ERN2/TMC5"),
    "5": ("Macrophages", "high", "CD163/MRC1/CTSB/FCER1G"),
    "6": ("AT1-like", "medium", "CAV1/HOPX/CAV2 with incomplete AGER/PDPN/AQP5 and retained SFTPB/LPCAT1/ABCA3"),
    "7": ("Endothelial cells", "high", "EPAS1/PECAM1/VWF/BTNL9"),
    "8": ("Endothelial cells", "high", "VWF/PTPRB/PECAM1/EPAS1"),
    "9": ("Macrophages", "high", "PPARG/MRC1/CD163/MSR1"),
    "10": ("Basal cells", "high", "EGFR/KRT15/COL7A1/TP63/KRT5"),
    "11": ("Smooth muscle / mural cells", "high", "MYH11/LMOD1/CARMN/PDGFRB/COL4A1"),
    "12": ("MT-high AT2-like", "medium", "MT genes with SFTPC/SFTPA2/SFTPB/ABCA3"),
    "13": ("T cells", "high", "PTPRC/CD247/BCL11B/ITK/DOCK2"),
    "14": ("Cycling cells", "medium", "DIAPH3/FANCI/MELK/RRM2/ECT2/ANLN; epithelial identity requires focused panel"),
    "15": ("Lymphatic endothelial cells", "high", "PROX1/FLT4/LYVE1/CCL21/MMRN1/RELN with pan-endothelial expression"),
    "16": ("NA", "low", "50-cell CYL-skewed unresolved epithelial-like cluster without a complete cycling signature"),
    "17": ("NA", "low", "20-cell ZCP-skewed unresolved stromal-like cluster with COL1A2/DCN/COL3A1"),
}

EPITHELIAL_FOCUS_CLUSTERS = ["1", "4", "6", "10", "12", "14", "16"]
EPITHELIAL_DOTPLOT_GROUPS = [
    "Secretory epithelial", "Secretory / mucous epithelial", "AT1-like", "Basal cells",
    "MT-high AT2-like", "Cycling cells", "NA",
]
RARE_REVIEW_CLUSTERS = ["15", "17"]
RARE_DOTPLOT_GROUPS = ["Lymphatic endothelial cells", "NA"]

if set(CELL_TYPE_ORDER) != {values[0] for values in CLUSTER_ANNOTATIONS.values()}:
    raise RuntimeError("DOTPLOT_MARKERS labels must exactly match the reviewed merged cell types")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cyl-zip", type=Path, default=DEFAULT_INPUTS["CYL"])
    parser.add_argument("--zcp-zip", type=Path, default=DEFAULT_INPUTS["ZCP"])
    parser.add_argument("--output-dir", type=Path, default=PROJECT_DIR / "Results/Scanpy/E_CYL_ZCP")
    parser.add_argument("--min-genes", type=int, default=200)
    parser.add_argument("--max-genes", type=int, default=6000)
    parser.add_argument("--min-counts", type=int, default=500)
    parser.add_argument("--max-mt-percent", type=float, default=5.0)
    parser.add_argument("--min-cells-per-gene", type=int, default=3)
    parser.add_argument("--doublet-rate-per-1000", type=float, default=0.004)
    parser.add_argument("--target-sum", type=float, default=1e4)
    parser.add_argument("--n-hvg", type=int, default=2000)
    parser.add_argument("--pca-components", type=int, default=50)
    parser.add_argument("--n-pcs", type=int, default=30)
    parser.add_argument("--n-neighbors", type=int, default=15)
    parser.add_argument("--leiden-resolution", type=float, default=0.8)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--allow-unreviewed-clusters", action="store_true")
    return parser.parse_args()


def log(message: str) -> None:
    print(message, flush=True)


def validate_args(args: argparse.Namespace) -> None:
    for path in (args.cyl_zip, args.zcp_zip):
        if not path.is_file():
            raise FileNotFoundError(path)
    if not (0 <= args.min_genes < args.max_genes):
        raise ValueError("Require 0 <= min_genes < max_genes")
    if args.min_counts < 0 or args.max_mt_percent < 0 or args.min_cells_per_gene < 1:
        raise ValueError("Invalid QC thresholds")
    if args.doublet_rate_per_1000 <= 0 or args.target_sum <= 0:
        raise ValueError("Doublet-rate coefficient and target sum must be positive")
    if min(args.n_hvg, args.pca_components, args.n_pcs, args.n_neighbors) < 1:
        raise ValueError("Dimensionality parameters must be positive")
    results_root = (PROJECT_DIR / "Results").resolve()
    output_dir = args.output_dir.resolve()
    if results_root not in output_dir.parents:
        raise ValueError(f"Output directory must be below {results_root}: {output_dir}")
    if output_dir.exists() and any(output_dir.iterdir()) and not args.overwrite:
        raise FileExistsError(f"Output directory is not empty: {output_dir}; use --overwrite to replace named files")


def prepare_output_dirs(output_dir: Path) -> dict[str, Path]:
    paths = {
        "root": output_dir,
        "tables": output_dir / "tables",
        "objects": output_dir / "objects",
        "qc_figures": output_dir / "figures/qc",
        "doublet_figures": output_dir / "figures/doublet",
        "dimension_figures": output_dir / "figures/dimensionality",
        "cluster_figures": output_dir / "figures/clustering",
        "annotation_figures": output_dir / "figures/annotation",
    }
    for path in paths.values():
        path.mkdir(parents=True, exist_ok=True)
    return paths


def save_figure(path: Path, figure=None) -> None:
    fig = figure if figure is not None else plt.gcf()
    fig.savefig(path, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def save_marker_dotplot(adata: ad.AnnData, marker_panels: dict[str, list[str]], groupby: str, path: Path) -> None:
    categories = list(adata.obs[groupby].cat.categories)
    if set(categories) != set(marker_panels):
        raise ValueError(f"Dotplot row categories and marker-panel labels differ: rows={categories}, panels={list(marker_panels)}")
    present_panels = {label: [gene for gene in marker_panels[label] if gene in adata.raw.var_names] for label in categories}
    empty_panels = [label for label, genes in present_panels.items() if not genes]
    if empty_panels:
        raise ValueError(f"No configured dotplot marker is present for: {empty_panels}")
    marker_order = [gene for genes in present_panels.values() for gene in genes]
    n_groups = int(adata.obs[groupby].nunique())
    figsize = (max(12, 0.32 * len(marker_order) + 5), max(5, 0.48 * n_groups + 3))
    sc.tl.dendrogram(
        adata,
        groupby=groupby,
        var_names=marker_order,
        use_raw=True,
        cor_method="pearson",
        linkage_method="complete",
    )
    dotplot = sc.pl.dotplot(
        adata,
        present_panels,
        groupby=groupby,
        use_raw=True,
        dendrogram=True,
        var_group_rotation=90,
        figsize=figsize,
        show=False,
        return_fig=True,
    )
    dotplot.style(
        cmap="Reds",
        dot_min=0,
        dot_max=0.6,
        smallest_dot=0,
        largest_dot=180,
        dot_edge_color="#888888",
        dot_edge_lw=0.5,
        size_exponent=1.5,
        grid=False,
    ).legend(
        size_title="Fraction of cells\nin group (%)",
        colorbar_title="Mean expression\nin group",
        width=1.8,
    )
    dotplot.savefig(path, dpi=220, bbox_inches="tight")
    plt.close("all")


def read_10x_zip(path: Path, cohort: str) -> ad.AnnData:
    with tempfile.TemporaryDirectory(prefix=f"{cohort}_10x_") as temporary_dir:
        with zipfile.ZipFile(path) as archive:
            archive.extractall(temporary_dir)
        matrix_dir = Path(temporary_dir) / "filtered_feature_bc_matrix"
        if not matrix_dir.is_dir():
            raise ValueError(f"Missing filtered_feature_bc_matrix in {path}")
        sample = sc.read_10x_mtx(matrix_dir, var_names="gene_symbols", make_unique=True)
    sample.obs_names = pd.Index([f"{cohort}_{barcode}" for barcode in sample.obs_names])
    sample.obs["cohort"] = cohort
    return sample


def add_qc_metrics(sample: ad.AnnData) -> None:
    sample.layers["counts"] = sample.X.copy()
    names = sample.var_names.str.upper()
    sample.var["mt"] = names.str.startswith("MT-")
    sample.var["ribo"] = names.str.startswith(("RPL", "RPS"))
    sc.pp.calculate_qc_metrics(sample, qc_vars=["mt", "ribo"], percent_top=None, log1p=False, inplace=True)


def plot_prefilter_qc(sample: ad.AnnData, cohort: str, figure_dir: Path) -> None:
    sc.pl.highest_expr_genes(sample, n_top=20, show=False)
    save_figure(figure_dir / f"{cohort}_highest_expression_genes.png")
    metrics = ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo"]
    fig, axes = plt.subplots(1, len(metrics), figsize=(16, 4), layout="constrained")
    for axis, metric in zip(axes, metrics):
        axis.violinplot(sample.obs[metric].dropna().to_numpy(), showmedians=True)
        axis.set_title(metric)
        axis.set_xticks([])
    fig.suptitle(f"{cohort}: pre-filter per-cell QC")
    save_figure(figure_dir / f"{cohort}_qc_violin_prefilter.png", fig)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), layout="constrained")
    sc.pl.scatter(sample, x="total_counts", y="pct_counts_mt", ax=axes[0], show=False)
    sc.pl.scatter(sample, x="total_counts", y="n_genes_by_counts", ax=axes[1], show=False)
    fig.suptitle(f"{cohort}: pre-filter QC relationships")
    save_figure(figure_dir / f"{cohort}_qc_scatter_prefilter.png", fig)


def rank_genes_table(adata: ad.AnnData) -> pd.DataFrame:
    try:
        return sc.get.rank_genes_groups_df(adata, group=None)
    except (AttributeError, TypeError):
        result = adata.uns["rank_genes_groups"]
        rows = []
        for group in result["names"].dtype.names:
            for rank, gene in enumerate(result["names"][group], start=1):
                row = {"group": group, "rank": rank, "names": gene}
                for key in ("scores", "logfoldchanges", "pvals", "pvals_adj"):
                    if key in result:
                        row[key] = result[key][group][rank - 1]
                rows.append(row)
        return pd.DataFrame(rows)


def write_json(path: Path, value: object) -> None:
    with path.open("w") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")


def main() -> None:
    args = parse_args()
    validate_args(args)
    paths = prepare_output_dirs(args.output_dir.resolve())
    sc.settings.verbosity = 2
    sc.set_figure_params(dpi=100, frameon=False)
    parameters = {
        "inputs": {"CYL": str(args.cyl_zip.resolve()), "ZCP": str(args.zcp_zip.resolve())},
        "min_genes": args.min_genes,
        "max_genes": args.max_genes,
        "min_counts": args.min_counts,
        "max_mt_percent": args.max_mt_percent,
        "min_cells_per_gene": args.min_cells_per_gene,
        "doublet_rate_formula": f"{args.doublet_rate_per_1000} * n_cells / 1000",
        "target_sum": args.target_sum,
        "n_hvg": args.n_hvg,
        "regress_covariates": False,
        "pca_components": args.pca_components,
        "n_pcs_neighbors": args.n_pcs,
        "n_neighbors": args.n_neighbors,
        "sample_key": "sample",
        "group_key": "group",
        "harmony_batch_key": "cohort",
        "harmony_basis": "X_pca_harmony",
        "before_umap_basis": "X_umap_before_harmony",
        "after_umap_basis": "X_umap_after_harmony",
        "harmony_max_iterations": 20,
        "leiden_resolution": args.leiden_resolution,
        "marker_test": {"groupby": "leiden", "method": "wilcoxon", "use_raw": True},
        "reviewed_cluster_ids": sorted(CLUSTER_ANNOTATIONS, key=int),
        "epithelial_focus_clusters": EPITHELIAL_FOCUS_CLUSTERS,
        "rare_review_clusters": RARE_REVIEW_CLUSTERS,
        "marker_panels": MARKER_GENES,
        "dotplot_markers": DOTPLOT_MARKERS,
        "dotplot_style": {
            "rows": "Merged reviewed cell type",
            "columns": "3-4 markers grouped by cell type",
            "shared_axis_order": CELL_TYPE_ORDER,
            "dendrogram": "Pearson/complete linkage on displayed raw log-normalized markers",
            "cmap": "Reds",
            "dot_min": 0,
            "dot_max": 0.6,
        },
        "seed": args.seed,
    }
    versions = {
        "python": platform.python_version(),
        "scanpy": sc.__version__,
        "anndata": ad.__version__,
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scrublet": getattr(scr, "__version__", "0.2.3"),
    }
    write_json(paths["root"] / "parameters.json", parameters)
    write_json(paths["root"] / "software_versions.json", versions)

    log("[1/8] Reading CYL and ZCP 10x matrices")
    samples = {
        "CYL": read_10x_zip(args.cyl_zip.resolve(), "CYL"),
        "ZCP": read_10x_zip(args.zcp_zip.resolve(), "ZCP"),
    }
    input_shapes = {cohort: [int(sample.n_obs), int(sample.n_vars)] for cohort, sample in samples.items()}

    log("[2/8] Calculating and exporting per-cohort QC")
    qc_descriptions = []
    prefilter_tables = []
    for cohort, sample in samples.items():
        add_qc_metrics(sample)
        sample.obs["pass_basic_qc"] = (
            (sample.obs["n_genes_by_counts"] >= args.min_genes)
            & (sample.obs["n_genes_by_counts"] < args.max_genes)
            & (sample.obs["total_counts"] >= args.min_counts)
            & (sample.obs["pct_counts_mt"] < args.max_mt_percent)
        )
        metrics = ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo"]
        description = sample.obs[metrics].describe(percentiles=[0.01, 0.05, 0.5, 0.95, 0.99])
        description.insert(0, "cohort", cohort)
        description.insert(1, "statistic", description.index)
        qc_descriptions.append(description.reset_index(drop=True))
        table = sample.obs[["cohort", *metrics, "pass_basic_qc"]].copy()
        table.index.name = "cell_id"
        prefilter_tables.append(table)
        plot_prefilter_qc(sample, cohort, paths["qc_figures"])
    pd.concat(qc_descriptions, ignore_index=True).to_csv(paths["tables"] / "qc_distribution_summary.tsv", sep="\t", index=False)
    pd.concat(prefilter_tables).to_csv(paths["tables"] / "cell_qc_prefilter.tsv.gz", sep="\t", compression="gzip")

    log("[3/8] Filtering each cohort and running Scrublet")
    retained_samples = {}
    doublet_tables = []
    qc_summary_rows = []
    for cohort, sample in samples.items():
        sample_qc = sample[sample.obs["pass_basic_qc"]].copy()
        sc.pp.filter_genes(sample_qc, min_cells=args.min_cells_per_gene)
        if sample_qc.n_obs < 100:
            raise RuntimeError(f"{cohort} has fewer than 100 cells after basic QC")
        expected_rate = args.doublet_rate_per_1000 * sample_qc.n_obs / 1000
        scrub = scr.Scrublet(sample_qc.layers["counts"], expected_doublet_rate=expected_rate, random_state=args.seed)
        scores, predictions = scrub.scrub_doublets(n_prin_comps=30, use_approx_neighbors=False, verbose=True)
        sample_qc.obs["doublet_score"] = scores
        sample_qc.obs["predicted_doublet"] = predictions
        scrub.plot_histogram()
        save_figure(paths["doublet_figures"] / f"{cohort}_scrublet_histogram.png")
        doublet_table = sample_qc.obs[
            ["cohort", "n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo", "doublet_score", "predicted_doublet"]
        ].copy()
        doublet_table.index.name = "cell_id"
        doublet_tables.append(doublet_table)
        retained = sample_qc[~sample_qc.obs["predicted_doublet"]].copy()
        retained_samples[cohort] = retained
        qc_summary_rows.append(
            {
                "cohort": cohort,
                "input_cells": sample.n_obs,
                "pass_basic_qc": sample_qc.n_obs,
                "expected_doublet_rate": expected_rate,
                "predicted_doublets": int(predictions.sum()),
                "retained_singlets": retained.n_obs,
            }
        )
    qc_summary = pd.DataFrame(qc_summary_rows).set_index("cohort")
    qc_summary.to_csv(paths["tables"] / "qc_doublet_summary.tsv", sep="\t")
    pd.concat(doublet_tables).to_csv(paths["tables"] / "scrublet_per_cell.tsv.gz", sep="\t", compression="gzip")

    log("[4/8] Merging singlets, normalizing, and selecting batch-aware HVGs")
    adata_qc = ad.concat(list(retained_samples.values()), join="outer", merge="same", index_unique=None)
    adata_qc.obs["cohort"] = adata_qc.obs["cohort"].astype("category")
    adata_qc.obs["sample"] = adata_qc.obs["cohort"].astype(str).astype("category")
    if adata_qc.obs_names.duplicated().any():
        raise RuntimeError("Duplicated cell IDs after cohort merge")
    sc.pp.normalize_total(adata_qc, target_sum=args.target_sum)
    sc.pp.log1p(adata_qc)
    adata_qc.raw = adata_qc
    sc.pp.highly_variable_genes(adata_qc, n_top_genes=args.n_hvg, batch_key="cohort", flavor="seurat")
    hvg_columns = [column for column in ["highly_variable", "highly_variable_nbatches", "means", "dispersions", "dispersions_norm"] if column in adata_qc.var]
    adata_qc.var[hvg_columns].to_csv(paths["tables"] / "highly_variable_genes.tsv", sep="\t", index_label="gene")
    sc.pl.highly_variable_genes(adata_qc, show=False)
    save_figure(paths["dimension_figures"] / "highly_variable_genes.png")

    log("[5/8] Running PCA and Harmony batch correction")
    adata_work = adata_qc[:, adata_qc.var["highly_variable"]].copy()
    sc.pp.scale(adata_work, max_value=10)
    n_pca = min(args.pca_components, adata_work.n_obs - 1, adata_work.n_vars - 1)
    if n_pca < args.n_pcs:
        raise RuntimeError(f"Only {n_pca} PCA components are possible, fewer than n_pcs={args.n_pcs}")
    sc.tl.pca(adata_work, n_comps=n_pca, svd_solver="arpack", random_state=args.seed)
    pd.DataFrame(
        {"pc": np.arange(1, len(adata_work.uns["pca"]["variance_ratio"]) + 1), "variance_ratio": adata_work.uns["pca"]["variance_ratio"]}
    ).to_csv(paths["tables"] / "pca_variance_ratio.tsv", sep="\t", index=False)
    sc.pl.pca_variance_ratio(adata_work, n_pcs=n_pca, log=True, show=False)
    save_figure(paths["dimension_figures"] / "pca_variance_ratio.png")
    sc.pl.pca(adata_work, color="cohort", show=False)
    save_figure(paths["dimension_figures"] / "pca_before_harmony_by_cohort.png")
    n_harmony_clusters = min(100, max(2, int(np.round(adata_work.n_obs / 30))))
    harmony_basis = "X_pca_harmony"
    sce.pp.harmony_integrate(
        adata_work,
        key="cohort",
        basis="X_pca",
        adjusted_basis=harmony_basis,
        nclust=n_harmony_clusters,
        sigma=np.repeat(0.1, n_harmony_clusters),
        max_iter_harmony=20,
        random_state=args.seed,
    )
    parameters["n_harmony_clusters"] = n_harmony_clusters

    log("[6/8] Building before/after-Harmony neighbours and UMAPs, then Leiden clusters")
    before_neighbors_key = "neighbors_before_harmony"
    before_umap_basis = "X_umap_before_harmony"
    after_umap_basis = "X_umap_after_harmony"
    sc.pp.neighbors(
        adata_work,
        n_neighbors=args.n_neighbors,
        use_rep="X_pca",
        n_pcs=args.n_pcs,
        random_state=args.seed,
        key_added=before_neighbors_key,
    )
    sc.tl.umap(adata_work, neighbors_key=before_neighbors_key, random_state=args.seed)
    adata_work.obsm[before_umap_basis] = adata_work.obsm["X_umap"].copy()
    sc.pp.neighbors(adata_work, n_neighbors=args.n_neighbors, use_rep=harmony_basis, n_pcs=args.n_pcs, random_state=args.seed)
    sc.tl.umap(adata_work, random_state=args.seed)
    adata_work.obsm[after_umap_basis] = adata_work.obsm["X_umap"].copy()
    sc.tl.leiden(adata_work, resolution=args.leiden_resolution, random_state=args.seed)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), layout="constrained")
    for axis, basis, title in zip(
        axes,
        ["umap_before_harmony", "umap_after_harmony"],
        ["Before Harmony", "After Harmony"],
    ):
        sc.pl.embedding(
            adata_work,
            basis=basis,
            color="sample",
            ax=axis,
            show=False,
            size=10,
            alpha=0.75,
            frameon=False,
            title=title,
            palette=["#0072B2", "#E69F00"],
        )
    save_figure(paths["cluster_figures"] / "umap_before_after_harmony_by_sample.png", fig)
    sc.pl.embedding(adata_work, basis="umap_after_harmony", color=["total_counts", "pct_counts_mt"], frameon=False, show=False)
    save_figure(paths["cluster_figures"] / "umap_harmony_qc.png")
    sc.pl.embedding(
        adata_work,
        basis="umap_after_harmony",
        color="leiden",
        legend_loc="on data",
        frameon=False,
        title="Leiden clusters",
        show=False,
    )
    save_figure(paths["cluster_figures"] / "umap_leiden.png")
    sc.pl.embedding(
        adata_work,
        basis="umap_after_harmony",
        color="sample",
        legend_loc="right margin",
        size=10,
        alpha=0.75,
        frameon=False,
        title="Sample",
        palette=["#0072B2", "#E69F00"],
        show=False,
    )
    save_figure(paths["cluster_figures"] / "umap_sample.png")
    if "group" in adata_work.obs and adata_work.obs["group"].notna().any():
        sc.pl.embedding(
            adata_work,
            basis="umap_after_harmony",
            color="group",
            legend_loc="right margin",
            size=10,
            alpha=0.75,
            frameon=False,
            title="Group",
            show=False,
        )
        save_figure(paths["cluster_figures"] / "umap_group.png")
        group_status = {"available": True, "column": "group", "figure": "figures/clustering/umap_group.png"}
    else:
        group_status = {"available": False, "column": "group", "reason": "No verified group metadata supplied"}
        log("Group UMAP skipped: no verified obs['group'] metadata")
    write_json(paths["root"] / "group_metadata_status.json", group_status)
    cluster_counts = pd.crosstab(adata_work.obs["leiden"], adata_work.obs["sample"])
    cluster_fractions = pd.crosstab(adata_work.obs["leiden"], adata_work.obs["sample"], normalize="index")
    cluster_counts.to_csv(paths["tables"] / "cluster_sample_counts.tsv", sep="\t")
    cluster_fractions.to_csv(paths["tables"] / "cluster_sample_fractions.tsv", sep="\t")
    cluster_qc = adata_work.obs.groupby("leiden", observed=True).agg(
        cells=("cohort", "size"),
        median_genes=("n_genes_by_counts", "median"),
        median_counts=("total_counts", "median"),
        median_pct_mt=("pct_counts_mt", "median"),
        median_pct_ribo=("pct_counts_ribo", "median"),
        median_doublet_score=("doublet_score", "median"),
    )
    cluster_qc.to_csv(paths["tables"] / "cluster_qc_summary.tsv", sep="\t")

    log("[7/8] Ranking markers and validating manual annotation")
    sc.tl.rank_genes_groups(adata_work, "leiden", method="wilcoxon", use_raw=True)
    rank_genes_table(adata_work).to_csv(paths["tables"] / "leiden_ranked_markers.tsv.gz", sep="\t", index=False, compression="gzip")
    pd.DataFrame(adata_work.uns["rank_genes_groups"]["names"]).head(50).to_csv(
        paths["tables"] / "leiden_top50_marker_names.tsv", sep="\t", index=False
    )
    sc.pl.rank_genes_groups(adata_work, n_genes=20, sharey=False, show=False)
    save_figure(paths["annotation_figures"] / "leiden_top20_ranked_markers.png")
    pd.DataFrame(
        [(module, gene, gene in adata_work.raw.var_names) for module, genes in MARKER_GENES.items() for gene in genes],
        columns=["module", "gene", "present"],
    ).to_csv(paths["tables"] / "annotation_marker_presence.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "cell_type_order": group_order,
                "cell_type": cell_type,
                "marker_order": marker_order,
                "gene": gene,
                "present": gene in adata_work.raw.var_names,
            }
            for group_order, (cell_type, genes) in enumerate(DOTPLOT_MARKERS.items(), start=1)
            for marker_order, gene in enumerate(genes, start=1)
        ]
    ).to_csv(paths["tables"] / "dotplot_marker_panel.tsv", sep="\t", index=False)
    observed_clusters = set(adata_work.obs["leiden"].astype(str).unique())
    reviewed_clusters = set(CLUSTER_ANNOTATIONS)
    missing_clusters = sorted(reviewed_clusters - observed_clusters, key=int)
    new_clusters = sorted(observed_clusters - reviewed_clusters, key=int)
    annotation_matches = not missing_clusters and not new_clusters
    cell_type_map = {cluster: values[0] for cluster, values in CLUSTER_ANNOTATIONS.items()}
    if annotation_matches:
        mapped_cell_types = adata_work.obs["leiden"].astype(str).map(cell_type_map)
    else:
        # Cluster IDs are run-specific. If the complete set changes, even a reused
        # numeric ID may represent different biology, so no reviewed label transfers.
        mapped_cell_types = pd.Series("Unassigned", index=adata_work.obs_names, dtype="object")
    epithelial_clusters = [cluster for cluster in EPITHELIAL_FOCUS_CLUSTERS if cluster in observed_clusters]
    rare_clusters = [cluster for cluster in RARE_REVIEW_CLUSTERS if cluster in observed_clusters]
    if annotation_matches:
        adata_work.obs["cell_type"] = pd.Categorical(mapped_cell_types, categories=CELL_TYPE_ORDER, ordered=True)
        save_marker_dotplot(
            adata_work,
            DOTPLOT_MARKERS,
            "cell_type",
            paths["annotation_figures"] / "annotation_marker_dotplot.png",
        )
        epithelial_subset = adata_work[adata_work.obs["leiden"].astype(str).isin(epithelial_clusters)].copy()
        epithelial_subset.obs["cell_type"] = epithelial_subset.obs["cell_type"].cat.remove_unused_categories()
        epithelial_markers = {name: DOTPLOT_MARKERS[name] for name in EPITHELIAL_DOTPLOT_GROUPS}
        save_marker_dotplot(
            epithelial_subset,
            epithelial_markers,
            "cell_type",
            paths["annotation_figures"] / "epithelial_focus_marker_dotplot.png",
        )
        rare_subset = adata_work[adata_work.obs["leiden"].astype(str).isin(rare_clusters)].copy()
        rare_subset.obs["cell_type"] = rare_subset.obs["cell_type"].cat.remove_unused_categories()
        rare_markers = {name: DOTPLOT_MARKERS[name] for name in RARE_DOTPLOT_GROUPS}
        save_marker_dotplot(
            rare_subset,
            rare_markers,
            "cell_type",
            paths["annotation_figures"] / "rare_cluster_marker_dotplot.png",
        )
    else:
        adata_work.obs["cell_type"] = mapped_cell_types.astype("category")
        log("Cell-type dotplots skipped because the observed Leiden IDs do not match the reviewed set")

    if rare_clusters:
        rare_qc = cluster_qc.loc[cluster_qc.index.astype(str).isin(rare_clusters)].copy()
        rare_qc = rare_qc.join(cluster_counts.add_suffix("_cells")).join(cluster_fractions.add_suffix("_fraction"))
        rare_qc.to_csv(paths["tables"] / "rare_cluster_qc_summary.tsv", sep="\t")

    annotation_table = pd.DataFrame(
        [
            {"leiden": cluster, "cell_type": values[0], "confidence": values[1], "marker_evidence": values[2]}
            for cluster, values in CLUSTER_ANNOTATIONS.items()
        ]
    ).set_index("leiden")
    annotation_table.to_csv(paths["tables"] / "reviewed_cluster_annotations.tsv", sep="\t")
    annotation_guard = {
        "status": "pass" if annotation_matches else "requires_review",
        "matches_reviewed_clusters": annotation_matches,
        "observed_clusters": sorted(observed_clusters, key=int),
        "reviewed_clusters": sorted(reviewed_clusters, key=int),
        "missing_reviewed_clusters": missing_clusters,
        "new_unreviewed_clusters": new_clusters,
    }
    write_json(paths["root"] / "annotation_guard_status.json", annotation_guard)
    validation_qc = cluster_qc.copy()
    validation_qc.index = validation_qc.index.astype(str)
    validation_counts = cluster_counts.copy().add_suffix("_cells")
    validation_counts.index = validation_counts.index.astype(str)
    validation_fractions = cluster_fractions.copy().add_suffix("_fraction")
    validation_fractions.index = validation_fractions.index.astype(str)
    validation_cluster_ids = sorted(reviewed_clusters | observed_clusters, key=int)
    validation_annotations = annotation_table.reindex(validation_cluster_ids).copy()
    validation_annotations["cell_type"] = validation_annotations["cell_type"].fillna("Unassigned")
    validation_annotations["confidence"] = validation_annotations["confidence"].fillna("unreviewed")
    validation_annotations["marker_evidence"] = validation_annotations["marker_evidence"].fillna(
        "New cluster in this run; inspect ranked markers, QC, sample composition, and competing lineage programs"
    )
    cluster_validation = validation_annotations.join(validation_qc).join(validation_counts).join(validation_fractions)
    cluster_validation.insert(3, "observed_in_run", cluster_validation.index.isin(observed_clusters))
    cluster_validation.to_csv(paths["tables"] / "cluster_annotation_validation.tsv", sep="\t")
    (
        annotation_table.reset_index()
        .groupby("cell_type", sort=False, observed=True)
        .agg(
            leiden_clusters=("leiden", lambda values: ",".join(sorted(values.astype(str), key=int))),
            cluster_count=("leiden", "size"),
            minimum_confidence=("confidence", lambda values: "low" if "low" in set(values) else ("medium" if "medium" in set(values) else "high")),
        )
        .reindex(CELL_TYPE_ORDER)
        .to_csv(paths["tables"] / "cell_type_cluster_membership.tsv", sep="\t")
    )
    if not annotation_matches:
        mismatch = {
            "message": "Leiden IDs differ from the reviewed 18-cluster notebook result; all cells are Unassigned.",
            "missing_reviewed_clusters": missing_clusters,
            "new_unreviewed_clusters": new_clusters,
        }
        write_json(paths["root"] / "annotation_mismatch.json", mismatch)
        log(f"WARNING: {mismatch}")
    if annotation_matches:
        sc.pl.embedding(
            adata_work,
            basis="umap_after_harmony",
            color="cell_type",
            legend_loc="right margin",
            size=10,
            alpha=0.85,
            frameon=False,
            title="Cell-type annotation",
            show=False,
        )
        save_figure(paths["annotation_figures"] / "umap_cell_type.png")
    else:
        log("Cell-type UMAP skipped because all labels remain Unassigned pending review")
    pd.crosstab(adata_work.obs["cell_type"], adata_work.obs["sample"]).to_csv(
        paths["tables"] / "cell_type_sample_counts.tsv", sep="\t"
    )
    pd.crosstab(adata_work.obs["cell_type"], adata_work.obs["sample"], normalize="columns").to_csv(
        paths["tables"] / "cell_type_sample_fractions.tsv", sep="\t"
    )

    log("[8/8] Saving annotated data and run metadata")
    adata_work.obsm["X_umap_harmony"] = adata_work.obsm[after_umap_basis].copy()
    adata_work.uns["scanpy_run_parameters"] = parameters
    adata_work.uns["cluster_annotation_evidence"] = annotation_table.reset_index().to_dict(orient="list")
    adata_work.uns["annotation_matches_reviewed_clusters"] = annotation_matches
    cell_metadata = adata_work.obs.copy()
    cell_metadata["umap_before_harmony_1"] = adata_work.obsm[before_umap_basis][:, 0]
    cell_metadata["umap_before_harmony_2"] = adata_work.obsm[before_umap_basis][:, 1]
    cell_metadata["umap_after_harmony_1"] = adata_work.obsm[after_umap_basis][:, 0]
    cell_metadata["umap_after_harmony_2"] = adata_work.obsm[after_umap_basis][:, 1]
    cell_metadata["umap_1"] = adata_work.obsm[after_umap_basis][:, 0]
    cell_metadata["umap_2"] = adata_work.obsm[after_umap_basis][:, 1]
    cell_metadata.index.name = "cell_id"
    cell_metadata.to_csv(paths["tables"] / "cell_metadata_qc_annotation.tsv.gz", sep="\t", compression="gzip")
    h5ad_name = "rna_e_cyl_zcp_annotated.h5ad" if annotation_matches else "rna_e_cyl_zcp_unreviewed_clusters.h5ad"
    h5ad_path = paths["objects"] / h5ad_name
    adata_work.write_h5ad(h5ad_path, compression="gzip")
    write_json(paths["root"] / "parameters.json", parameters)
    summary = {
        "status": "complete" if annotation_matches else "annotation_requires_review",
        "input_shapes": input_shapes,
        "qc": qc_summary.reset_index().to_dict(orient="records"),
        "merged_singlets": int(adata_work.n_obs),
        "full_genes_in_raw": int(adata_work.raw.n_vars),
        "highly_variable_genes": int(adata_work.n_vars),
        "leiden_clusters": int(adata_work.obs["leiden"].nunique()),
        "group_metadata_available": bool(group_status["available"]),
        "annotation_matches_reviewed_clusters": annotation_matches,
        "new_unreviewed_clusters": new_clusters,
        "missing_reviewed_clusters": missing_clusters,
        "h5ad": str(h5ad_path),
    }
    write_json(paths["root"] / "run_summary.json", summary)
    log(json.dumps(summary, sort_keys=True))
    if not annotation_matches and not args.allow_unreviewed_clusters:
        raise RuntimeError(
            "Outputs were saved, but cluster IDs differ from the reviewed notebook annotation. "
            "Review ranked markers, cluster QC/sample composition, and lineage-marker evidence before updating CLUSTER_ANNOTATIONS."
        )
    log(f"Scanpy analysis complete: {h5ad_path}")


if __name__ == "__main__":
    main()
