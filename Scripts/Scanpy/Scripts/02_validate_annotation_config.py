#!/usr/bin/env python3
"""Validate the reviewed Scanpy cluster annotation and dotplot configuration."""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path


DEFAULT_PRODUCTION_SCRIPT = Path(__file__).with_name("01_run_e_scanpy.py")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--production-script", type=Path, default=DEFAULT_PRODUCTION_SCRIPT)
    parser.add_argument("--json-out", type=Path, help="Optional path for the same machine-readable report printed to stdout")
    return parser.parse_args()


def load_production_module(path: Path):
    if not path.is_file():
        raise FileNotFoundError(path)
    spec = importlib.util.spec_from_file_location("ipf_scanpy_production", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load production script: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def validate(module) -> dict[str, object]:
    errors: list[str] = []
    annotations = module.CLUSTER_ANNOTATIONS
    cluster_ids = sorted(annotations, key=int)
    expected_contiguous_ids = [str(value) for value in range(len(cluster_ids))]
    if cluster_ids != expected_contiguous_ids:
        errors.append(f"Reviewed cluster IDs are not contiguous from 0: {cluster_ids}")

    invalid_confidence = {
        cluster: values[1]
        for cluster, values in annotations.items()
        if values[1] not in {"high", "medium", "low"}
    }
    if invalid_confidence:
        errors.append(f"Invalid confidence values: {invalid_confidence}")

    annotation_cell_types = {values[0] for values in annotations.values()}
    marker_cell_types = set(module.DOTPLOT_MARKERS)
    if annotation_cell_types != marker_cell_types:
        errors.append(
            "Reviewed cell types and dotplot groups differ: "
            f"annotations={sorted(annotation_cell_types)}, markers={sorted(marker_cell_types)}"
        )
    if list(module.DOTPLOT_MARKERS) != list(module.CELL_TYPE_ORDER):
        errors.append("CELL_TYPE_ORDER must be identical to DOTPLOT_MARKERS insertion order")

    marker_count_errors = {
        cell_type: len(genes)
        for cell_type, genes in module.DOTPLOT_MARKERS.items()
        if not 3 <= len(genes) <= 4
    }
    if marker_count_errors:
        errors.append(f"Dotplot groups must contain 3-4 markers: {marker_count_errors}")

    duplicate_markers: dict[str, list[str]] = {}
    marker_owners: dict[str, list[str]] = {}
    for cell_type, genes in module.DOTPLOT_MARKERS.items():
        if len(genes) != len(set(genes)):
            errors.append(f"Duplicate marker inside {cell_type}: {genes}")
        for gene in genes:
            marker_owners.setdefault(gene, []).append(cell_type)
    duplicate_markers = {gene: owners for gene, owners in marker_owners.items() if len(owners) > 1}
    if duplicate_markers:
        errors.append(f"Markers assigned to multiple dotplot groups: {duplicate_markers}")

    reviewed_set = set(cluster_ids)
    focus_set = set(module.EPITHELIAL_FOCUS_CLUSTERS)
    rare_set = set(module.RARE_REVIEW_CLUSTERS)
    if not focus_set <= reviewed_set:
        errors.append(f"Unknown epithelial-focus clusters: {sorted(focus_set - reviewed_set, key=int)}")
    if not rare_set <= reviewed_set:
        errors.append(f"Unknown rare-review clusters: {sorted(rare_set - reviewed_set, key=int)}")

    focus_cell_types = {annotations[cluster][0] for cluster in module.EPITHELIAL_FOCUS_CLUSTERS}
    if focus_cell_types != set(module.EPITHELIAL_DOTPLOT_GROUPS):
        errors.append(
            "Epithelial focus cluster labels and dotplot groups differ: "
            f"clusters={sorted(focus_cell_types)}, groups={sorted(module.EPITHELIAL_DOTPLOT_GROUPS)}"
        )
    rare_cell_types = {annotations[cluster][0] for cluster in module.RARE_REVIEW_CLUSTERS}
    if rare_cell_types != set(module.RARE_DOTPLOT_GROUPS):
        errors.append(
            "Rare-review cluster labels and dotplot groups differ: "
            f"clusters={sorted(rare_cell_types)}, groups={sorted(module.RARE_DOTPLOT_GROUPS)}"
        )

    merged_membership: dict[str, list[str]] = {}
    for cluster in cluster_ids:
        merged_membership.setdefault(annotations[cluster][0], []).append(cluster)

    return {
        "status": "pass" if not errors else "fail",
        "production_script": str(Path(module.__file__).resolve()),
        "reviewed_cluster_ids": cluster_ids,
        "reviewed_cluster_count": len(cluster_ids),
        "merged_cell_type_count": len(merged_membership),
        "merged_cell_type_membership": merged_membership,
        "epithelial_focus_clusters": list(module.EPITHELIAL_FOCUS_CLUSTERS),
        "rare_review_clusters": list(module.RARE_REVIEW_CLUSTERS),
        "dotplot_marker_counts": {cell_type: len(genes) for cell_type, genes in module.DOTPLOT_MARKERS.items()},
        "errors": errors,
    }


def main() -> None:
    args = parse_args()
    module = load_production_module(args.production_script.resolve())
    report = validate(module)
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    print(payload, end="")
    if args.json_out is not None:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(payload)
    if report["status"] != "pass":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
