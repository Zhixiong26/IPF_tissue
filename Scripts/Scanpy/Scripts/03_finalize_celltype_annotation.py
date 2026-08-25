#!/usr/bin/env python3
"""Finalize the reviewed 2026-08-25 CYL/ZCP cell-type annotation.

The production run yielded Leiden clusters 0-16.  Cluster 9 contains several
biologically distinct populations, so this script deterministically reclusters
only that cluster before assigning reviewed labels.  It never overwrites the
unreviewed source object.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc


PROJECT_DIR = Path("/home/lijia/luozhixiong/IPF_tissue")
DEFAULT_INPUT = (
    PROJECT_DIR
    / "Results/Scanpy/E_CYL_ZCP_20260825/objects/rna_e_cyl_zcp_unreviewed_clusters.h5ad"
)
DEFAULT_OUTPUT = (
    PROJECT_DIR
    / "Results/Scanpy/E_CYL_ZCP_20260825/reviewed_annotation_20260825"
)

# Labels for clusters that do not require local reclustering.
CLUSTER_MAP = {
    "0": ("AT2", "high", "SFTPC/SFTPB/ABCA3/LPCAT1"),
    "1": ("Secretory epithelial", "medium", "NEDD4L/SFTPB/SFTA3"),
    "2": ("Fibroblasts", "high", "COL1A2/COL3A1/COL1A1/DCN"),
    "3": ("Ciliated cells", "high", "CFAP/DNAH/HYDIN"),
    "4": ("Macrophages", "high", "CD163/MRC1/PPARG/DOCK2"),
    "5": ("AT1-like", "medium", "CAV1/HOPX with retained surfactant genes"),
    "6": ("Secretory / mucous epithelial", "high", "BPIFB1/MUC4/ERN2/TMC5"),
    "7": ("Endothelial cells", "high", "EPAS1/PECAM1/VWF"),
    "8": ("Endothelial cells", "high", "VWF/PTPRB/PECAM1"),
    "10": ("Basal cells", "high", "EGFR/KRT15/TP63/KRT5"),
    "11": ("Smooth muscle / mural cells", "high", "MYH11/LMOD1/CARMN/PDGFRB"),
    "12": ("MT-high AT2-like", "medium", "MT genes with SFTPC/SFTPA2/SFTPB"),
    "13": ("Cycling cells", "high", "DIAPH3/FANCI/RRM2/ANLN/TOP2A/MKI67"),
    "14": (
        "Lymphatic endothelial cells",
        "high",
        "MMRN1/RELN/FLT4/PROX1/LYVE1/CCL21",
    ),
    "15": ("NA", "low", "50-cell CYL-skewed unresolved epithelial-like population"),
    "16": ("NA", "low", "20-cell ZCP-skewed COL11A1/CEMIP stromal-like population"),
}

# Subcluster identifiers are stable for the fixed input, seed, PCA and parameters.
CLUSTER9_MAP = {
    "0": ("Secretory epithelial", "medium", "SFTPB/SFTA3/SCNN1A/TMC5"),
    "1": ("Macrophages", "high", "CD163/MRC1/MERTK/F13A1"),
    "2": ("T cells", "high", "PTPRC/ITK/CD247/BCL11B/SKAP1/IL7R"),
    "3": ("Plasma cells", "medium", "IGHG1/IGHG3/XBP1/FCRL5"),
    "4": ("Mast cells", "high", "KIT/CPA3/HDC/MS4A2/TPSB2"),
}

SUBCLUSTER_MARKER_GUARD = {
    "0": {"SFTPB", "SFTA3", "SCNN1A", "TMC5"},
    "1": {"CD163", "MRC1", "MERTK", "F13A1"},
    "2": {"PTPRC", "ITK", "CD247", "BCL11B", "SKAP1", "IL7R"},
    "3": {"IGHG1", "IGHG3", "XBP1", "FCRL5"},
    "4": {"KIT", "CPA3", "HDC", "MS4A2", "TPSB2"},
}

CELL_TYPE_ORDER = [
    "AT2",
    "AT1-like",
    "Secretory epithelial",
    "Secretory / mucous epithelial",
    "Ciliated cells",
    "Basal cells",
    "Fibroblasts",
    "Smooth muscle / mural cells",
    "Endothelial cells",
    "Lymphatic endothelial cells",
    "Macrophages",
    "T cells",
    "Plasma cells",
    "Mast cells",
    "Cycling cells",
    "MT-high AT2-like",
    "NA",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-h5ad", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--n-neighbors", type=int, default=15)
    parser.add_argument("--resolution", type=float, default=0.25)
    parser.add_argument(
        "--plot-only",
        action="store_true",
        help="Plot an already-created reviewed H5AD in OUTPUT_DIR without rerunning annotation.",
    )
    parser.add_argument(
        "--publish-figure",
        type=Path,
        help="Optional second PNG path, for example the parent run's figures/annotation directory.",
    )
    parser.add_argument("--overwrite-figure", action="store_true")
    return parser.parse_args()


def annotation_review_notice(input_h5ad: Path) -> dict:
    guard_path = input_h5ad.parent.parent / "annotation_guard_status.json"
    notice = {
        "source_guard_json": str(guard_path),
        "source_guard_found": guard_path.is_file(),
        "source_annotation_guard_mismatch": False,
        "warning": None,
    }
    if guard_path.is_file():
        guard = json.loads(guard_path.read_text())
        notice["source_annotation_guard_mismatch"] = not bool(
            guard.get("matches_reviewed_clusters", False)
        )
        notice["source_guard_status"] = guard.get("status")
    if notice["source_annotation_guard_mismatch"]:
        notice["warning"] = (
            "The original run failed its automatic cluster-ID annotation guard. "
            "This UMAP is a reviewed post-run derivative, not an automatically transferred annotation."
        )
        print(f"WARNING: {notice['warning']}", flush=True)
    return notice


def plot_reviewed_umap(
    adata,
    output_png: Path,
    publish_png: Path | None,
    overwrite: bool,
) -> list[str]:
    if "cell_type" not in adata.obs:
        raise KeyError("Reviewed H5AD lacks obs['cell_type']")
    if "X_umap_after_harmony" not in adata.obsm:
        raise KeyError("Reviewed H5AD lacks obsm['X_umap_after_harmony']")
    labels = adata.obs["cell_type"].astype(str)
    if labels.eq("Unassigned").any() or labels.isna().any():
        raise RuntimeError("Cell-type UMAP requires reviewed labels; Unassigned/missing labels remain")
    observed_types = set(labels)
    unexpected_types = observed_types - set(CELL_TYPE_ORDER)
    if unexpected_types:
        raise RuntimeError(f"Unexpected reviewed cell types: {sorted(unexpected_types)}")

    present_order = [label for label in CELL_TYPE_ORDER if label in observed_types]
    adata.obs["cell_type"] = pd.Categorical(labels, categories=present_order, ordered=True)
    tab20 = plt.get_cmap("tab20")
    palette = ["#B0B0B0" if label == "NA" else tab20(i % 20) for i, label in enumerate(present_order)]

    targets = [output_png]
    if publish_png is not None and publish_png != output_png:
        targets.append(publish_png)
    for target in targets:
        if target.exists() and not overwrite:
            raise FileExistsError(f"Refusing to overwrite existing figure: {target}")
        target.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(12, 8), layout="constrained")
    sc.pl.embedding(
        adata,
        basis="umap_after_harmony",
        color="cell_type",
        palette=palette,
        legend_loc="right margin",
        legend_fontsize=8,
        size=12,
        alpha=0.85,
        frameon=False,
        title="CYL + ZCP reviewed cell types",
        ax=ax,
        show=False,
    )
    for target in targets:
        fig.savefig(target, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return [str(target) for target in targets]


def main() -> None:
    args = parse_args()
    if not args.input_h5ad.is_file():
        raise FileNotFoundError(args.input_h5ad)
    output_h5ad = args.output_dir / "objects/rna_e_cyl_zcp_annotated_reviewed.h5ad"
    canonical_figure = args.output_dir / "figures/annotation/umap_cell_type.png"
    notice = annotation_review_notice(args.input_h5ad)
    if args.plot_only:
        if not output_h5ad.is_file():
            raise FileNotFoundError(output_h5ad)
        adata = sc.read_h5ad(output_h5ad)
        figure_paths = plot_reviewed_umap(
            adata,
            canonical_figure,
            args.publish_figure,
            args.overwrite_figure,
        )
        status = {
            "status": "reviewed_cell_type_umap_complete",
            "reviewed_h5ad": str(output_h5ad),
            "n_cells": int(adata.n_obs),
            "n_cell_types": int(adata.obs["cell_type"].nunique()),
            "figure_paths": figure_paths,
            "annotation_review_notice": notice,
        }
        status_path = args.output_dir / "umap_cell_type_plot_status.json"
        status_path.write_text(json.dumps(status, indent=2, ensure_ascii=False) + "\n")
        print(json.dumps(status, indent=2, ensure_ascii=False), flush=True)
        return
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise FileExistsError(f"Refusing non-empty output directory: {args.output_dir}")

    tables = args.output_dir / "tables"
    objects = args.output_dir / "objects"
    tables.mkdir(parents=True, exist_ok=True)
    objects.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(args.input_h5ad)
    observed = set(adata.obs["leiden"].astype(str))
    expected = set(CLUSTER_MAP) | {"9"}
    if observed != expected:
        raise RuntimeError(
            f"Leiden guard failed: observed={sorted(observed, key=int)}, "
            f"expected={sorted(expected, key=int)}"
        )

    cluster9_mask = adata.obs["leiden"].astype(str).eq("9")
    cluster9 = adata[cluster9_mask].copy()
    sc.pp.neighbors(
        cluster9,
        n_neighbors=args.n_neighbors,
        use_rep="X_pca",
        random_state=args.seed,
    )
    sc.tl.leiden(
        cluster9,
        resolution=args.resolution,
        key_added="cluster9_subcluster",
        random_state=args.seed,
    )
    found_subclusters = set(cluster9.obs["cluster9_subcluster"].astype(str))
    if found_subclusters != set(CLUSTER9_MAP):
        raise RuntimeError(
            f"Cluster-9 guard failed: observed={sorted(found_subclusters)}, "
            f"expected={sorted(CLUSTER9_MAP)}"
        )

    cluster9.uns.setdefault("log1p", {})["base"] = None
    sc.tl.rank_genes_groups(
        cluster9,
        groupby="cluster9_subcluster",
        method="wilcoxon",
        pts=True,
        use_raw=True,
    )
    ranked = sc.get.rank_genes_groups_df(cluster9, group=None)
    for group, required in SUBCLUSTER_MARKER_GUARD.items():
        top100 = set(ranked.loc[ranked["group"].astype(str).eq(group), "names"].head(100))
        hits = required & top100
        if len(hits) < 2:
            raise RuntimeError(
                f"Marker guard failed for cluster 9.{group}: "
                f"only found {sorted(hits)} among expected {sorted(required)}"
            )
    ranked.to_csv(tables / "cluster9_subcluster_ranked_markers.tsv.gz", sep="\t", index=False)

    annotation_cluster = adata.obs["leiden"].astype(str).copy()
    annotation_cluster.loc[cluster9.obs_names] = (
        "9." + cluster9.obs["cluster9_subcluster"].astype(str)
    )
    adata.obs["annotation_cluster"] = annotation_cluster

    reviewed_map = dict(CLUSTER_MAP)
    reviewed_map.update({f"9.{key}": value for key, value in CLUSTER9_MAP.items()})
    for field, position in (("cell_type", 0), ("annotation_confidence", 1), ("marker_evidence", 2)):
        adata.obs[field] = adata.obs["annotation_cluster"].map(
            {key: value[position] for key, value in reviewed_map.items()}
        )
        if adata.obs[field].isna().any():
            raise RuntimeError(f"Unmapped cells remain in {field}")

    mapping_rows = [
        {
            "annotation_cluster": key,
            "source_leiden": key.split(".")[0],
            "cell_type": value[0],
            "confidence": value[1],
            "marker_evidence": value[2],
            "n_cells": int(adata.obs["annotation_cluster"].eq(key).sum()),
        }
        for key, value in reviewed_map.items()
    ]
    mapping = pd.DataFrame(mapping_rows).sort_values(
        "annotation_cluster", key=lambda x: x.str.split(".").map(lambda y: (int(y[0]), int(y[1]) if len(y) > 1 else -1))
    )
    mapping.to_csv(tables / "annotation_cluster_map.tsv", sep="\t", index=False)

    metadata_columns = [
        column
        for column in (
            "cohort", "sample", "leiden", "annotation_cluster", "cell_type",
            "annotation_confidence", "marker_evidence", "n_genes_by_counts",
            "total_counts", "pct_counts_mt", "doublet_score",
        )
        if column in adata.obs
    ]
    metadata = adata.obs[metadata_columns].copy()
    metadata.index.name = "cell_id"
    metadata.to_csv(tables / "celltype_annotations.tsv.gz", sep="\t")
    pd.crosstab(adata.obs["cell_type"], adata.obs["sample"]).to_csv(
        tables / "cell_type_sample_counts.tsv", sep="\t"
    )

    figure_paths = plot_reviewed_umap(
        adata,
        canonical_figure,
        args.publish_figure,
        args.overwrite_figure,
    )
    adata.write_h5ad(output_h5ad, compression="gzip")
    summary = {
        "status": "reviewed_annotation_complete",
        "source_h5ad": str(args.input_h5ad),
        "output_h5ad": str(output_h5ad),
        "celltype_table": str(tables / "celltype_annotations.tsv.gz"),
        "cell_type_umap": figure_paths,
        "annotation_review_notice": notice,
        "n_cells": int(adata.n_obs),
        "n_source_leiden_clusters": int(adata.obs["leiden"].nunique()),
        "n_annotation_clusters": int(adata.obs["annotation_cluster"].nunique()),
        "n_cell_types": int(adata.obs["cell_type"].nunique()),
        "cluster9_subcluster_counts": {
            str(key): int(value)
            for key, value in cluster9.obs["cluster9_subcluster"].value_counts().sort_index().items()
        },
        "parameters": {
            "seed": args.seed,
            "n_neighbors": args.n_neighbors,
            "resolution": args.resolution,
        },
    }
    (args.output_dir / "run_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
