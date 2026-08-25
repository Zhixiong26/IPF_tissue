#!/usr/bin/env python3
"""Select annotated cells and convert six-column Bismark coverage to ALLC."""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import hashlib
import json
import os
import subprocess
from pathlib import Path

import pandas as pd


COV_SUFFIX = "_allc.gz.cov"


def cell_id(path: Path) -> str:
    if not path.name.endswith(COV_SUFFIX):
        raise ValueError(f"Unexpected coverage filename: {path.name}")
    return path.name[: -len(COV_SUFFIX)]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def valid_allc(path: Path) -> bool:
    index = Path(f"{path}.tbi")
    return (
        path.is_file() and path.stat().st_size > 0 and os.access(path, os.R_OK)
        and index.is_file() and index.stat().st_size > 0 and os.access(index, os.R_OK)
    )


def convert_one(task: tuple[str, str, str, str]) -> dict[str, object]:
    cov_string, output_string, bgzip, tabix = task
    cov, output = Path(cov_string), Path(output_string)
    if valid_allc(output):
        return {"cell_id": cell_id(cov), "status": "reused"}

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.tmp.{os.getpid()}")
    temporary_index = Path(f"{temporary}.tbi")
    rows = 0
    previous: tuple[str, int] | None = None
    process = None
    try:
        with temporary.open("wb") as compressed:
            process = subprocess.Popen(
                [bgzip, "-@", "1", "-c"], stdin=subprocess.PIPE,
                stdout=compressed, stderr=subprocess.PIPE,
            )
            assert process.stdin is not None
            with cov.open("rt") as source:
                for line_number, line in enumerate(source, start=1):
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) != 6:
                        raise ValueError(f"{cov}:{line_number}: expected 6 columns")
                    chrom, start_s, end_s, _ratio, mc_s, uc_s = fields
                    try:
                        start, end, mc, uc = int(start_s), int(end_s), int(mc_s), int(uc_s)
                    except ValueError as exc:
                        raise ValueError(f"{cov}:{line_number}: non-integer coordinate/count") from exc
                    if start < 1 or end < start or mc < 0 or uc < 0:
                        raise ValueError(f"{cov}:{line_number}: invalid coordinate/count")
                    current = (chrom, start)
                    if current == previous:
                        raise ValueError(f"{cov}:{line_number}: duplicate CpG coordinate {chrom}:{start}")
                    previous = current
                    process.stdin.write(
                        f"{chrom}\t{start}\t+\tCGN\t{mc}\t{mc + uc}\t1\n".encode()
                    )
                    rows += 1
            process.stdin.close()
            stderr = process.stderr.read().decode(errors="replace") if process.stderr else ""
            if process.wait():
                raise RuntimeError(f"bgzip failed for {cov}: {stderr.strip()}")
        subprocess.run(
            [tabix, "-f", "-s", "1", "-b", "2", "-e", "2", str(temporary)],
            check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True,
        )
        temporary.replace(output)
        temporary_index.replace(Path(f"{output}.tbi"))
    except Exception:
        if process is not None and process.poll() is None:
            process.kill()
        temporary.unlink(missing_ok=True)
        temporary_index.unlink(missing_ok=True)
        raise
    return {"cell_id": cell_id(cov), "status": "built", "sites": rows}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cov-dir", type=Path, default=Path(os.environ["IPF_COV_DIR"]))
    parser.add_argument("--annotation", type=Path, default=Path(os.environ["IPF_ANNOTATION"]))
    parser.add_argument("--output-dir", type=Path, default=Path(os.environ["IPF_ALLC_DIR"]))
    parser.add_argument(
        "--existing-allc-dir", type=Path,
        default=Path(os.environ["IPF_EXISTING_ALLC_DIR"]),
        help="Reuse readable cell.allc.tsv.gz files from this directory before converting cov",
    )
    parser.add_argument("--allc-table", type=Path, default=Path(os.environ["IPF_ALLC_TABLE"]))
    parser.add_argument("--threads", type=int, default=int(os.environ["IPF_THREADS"]))
    parser.add_argument("--expected-cells", type=int, default=int(os.environ["IPF_EXPECTED_CELLS"]))
    parser.add_argument("--max-cells", type=int, default=int(os.environ.get("IPF_MAX_CELLS", "0")))
    parser.add_argument(
        "--balanced-cohorts", action="store_true",
        default=os.environ.get("IPF_BALANCED_COHORTS") == "1",
        help="For max-cells tests, select cells round-robin across filename prefixes",
    )
    parser.add_argument(
        "--include-unannotated", action="store_true",
        default=os.environ.get("IPF_INCLUDE_UNANNOTATED") == "1",
    )
    parser.add_argument("--verify-only", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.threads < 1 or args.expected_cells < 1 or args.max_cells < 0:
        raise ValueError("threads/expected-cells must be positive and max-cells non-negative")
    files = sorted(args.cov_dir.glob(f"*{COV_SUFFIX}"))
    if not files:
        raise FileNotFoundError(f"No *{COV_SUFFIX} files in {args.cov_dir}")
    annotation = pd.read_csv(args.annotation, sep="\t", dtype=str)
    required = {"cell_id", "manual_celltype"}
    if not required.issubset(annotation.columns):
        raise ValueError(f"Annotation needs columns {sorted(required)}")
    if annotation["cell_id"].duplicated().any():
        raise ValueError("Annotation contains duplicate cell_id values")
    annotation = annotation.set_index("cell_id")
    selected: list[tuple[str, Path]] = []
    for path in files:
        name = cell_id(path)
        annotated = name in annotation.index
        if args.include_unannotated or annotated:
            selected.append((name, path.resolve()))
    if args.max_cells:
        if args.balanced_cohorts:
            groups: dict[str, list[tuple[str, Path]]] = {}
            for item in selected:
                groups.setdefault(item[0].split("_", 1)[0], []).append(item)
            balanced: list[tuple[str, Path]] = []
            while len(balanced) < args.max_cells and any(groups.values()):
                for cohort in sorted(groups):
                    if groups[cohort] and len(balanced) < args.max_cells:
                        balanced.append(groups[cohort].pop(0))
            selected = balanced
        else:
            selected = selected[: args.max_cells]
    if not selected:
        raise RuntimeError("No cells remain after annotation filtering")
    if not args.max_cells and len(selected) != args.expected_cells:
        raise RuntimeError(f"Expected {args.expected_cells:,} selected cells, found {len(selected):,}")

    summary = {
        "coverage_files": len(files), "selected_cells": len(selected),
        "excluded_cells": len(files) - len(selected),
        "annotation_sha256": sha256(args.annotation),
        "cell_ids_matched_to_annotation": sum(name in annotation.index for name, _path in selected),
        "include_unannotated": args.include_unannotated, "max_cells": args.max_cells,
        "balanced_cohorts": args.balanced_cohorts,
        "cohorts": pd.Series([name.split("_", 1)[0] for name, _path in selected]).value_counts().to_dict(),
    }
    summary["existing_allc_available"] = sum(
        valid_allc(args.existing_allc_dir / f"{name}.allc.tsv.gz")
        for name, _path in selected
    )
    print(json.dumps(summary, indent=2), flush=True)
    if args.verify_only:
        return

    env_prefix = Path(os.environ["IPF_ALLCOOLS_ENV"])
    bgzip, tabix = str(env_prefix / "bin/bgzip"), str(env_prefix / "bin/tabix")
    for executable in (bgzip, tabix):
        if not os.access(executable, os.X_OK):
            raise FileNotFoundError(executable)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.allc_table.parent.mkdir(parents=True, exist_ok=True)
    source_manifest = args.allc_table.parent / "source_cov_manifest.tsv"
    source_lines = [
        f"{name}\t{path}\t{path.stat().st_size}\t{path.stat().st_mtime_ns}\n"
        for name, path in selected
    ]
    source_text = "cell_id\tpath\tsize\tmtime_ns\n" + "".join(source_lines)
    if source_manifest.exists() and source_manifest.read_text() != source_text:
        raise RuntimeError(
            f"Coverage source manifest changed; use a new output root: {source_manifest.parent}"
        )
    if not source_manifest.exists() and any(args.output_dir.glob("*.allc.tsv.gz")):
        raise RuntimeError(f"Unversioned ALLC files found in {args.output_dir}")
    if not source_manifest.exists():
        temporary_manifest = source_manifest.with_suffix(".tmp.tsv")
        temporary_manifest.write_text(source_text)
        temporary_manifest.replace(source_manifest)
    allc_paths: dict[str, Path] = {}
    tasks = []
    reused_local = reused_existing = 0
    for name, path in selected:
        local_allc = args.output_dir / f"{name}.allc.tsv.gz"
        existing_allc = args.existing_allc_dir / f"{name}.allc.tsv.gz"
        if valid_allc(local_allc):
            allc_paths[name] = local_allc.resolve()
            reused_local += 1
        elif valid_allc(existing_allc):
            allc_paths[name] = existing_allc.resolve()
            reused_existing += 1
        else:
            tasks.append((str(path), str(local_allc), bgzip, tabix))

    built = 0
    if tasks:
        with cf.ProcessPoolExecutor(max_workers=args.threads) as executor:
            for completed, result in enumerate(executor.map(convert_one, tasks), start=1):
                built += result["status"] == "built"
                name = str(result["cell_id"])
                allc_paths[name] = (args.output_dir / f"{name}.allc.tsv.gz").resolve()
                if completed % 25 == 0 or completed == len(tasks):
                    print(
                        f"ALLC conversion {completed:,}/{len(tasks):,} (built={built:,})",
                        flush=True,
                    )
    print(
        f"ALLC sources: built={built:,}, local={reused_local:,}, "
        f"existing={reused_existing:,}", flush=True,
    )

    table_lines = []
    for name, _path in selected:
        allc = allc_paths.get(name)
        if allc is None:
            raise RuntimeError(f"No ALLC source selected for {name}")
        if not valid_allc(allc):
            raise RuntimeError(f"Missing/invalid ALLC: {allc}")
        table_lines.append(f"{name}\t{allc}\n")

    allc_manifest = args.allc_table.parent / "input_allc_manifest.tsv"
    allc_manifest_text = "cell_id\tpath\tsize\tmtime_ns\tindex_size\tindex_mtime_ns\n" + "".join(
        f"{name}\t{allc_paths[name]}\t{allc_paths[name].stat().st_size}\t"
        f"{allc_paths[name].stat().st_mtime_ns}\t"
        f"{Path(f'{allc_paths[name]}.tbi').stat().st_size}\t"
        f"{Path(f'{allc_paths[name]}.tbi').stat().st_mtime_ns}\n"
        for name, _path in selected
    )
    if allc_manifest.exists() and allc_manifest.read_text() != allc_manifest_text:
        raise RuntimeError(
            f"ALLC input manifest changed; use a new output root: {allc_manifest.parent}"
        )
    if not allc_manifest.exists():
        temporary_allc_manifest = allc_manifest.with_suffix(".tmp.tsv")
        temporary_allc_manifest.write_text(allc_manifest_text)
        temporary_allc_manifest.replace(allc_manifest)

    temporary = args.allc_table.with_suffix(".tmp.tsv")
    temporary.write_text("".join(table_lines))
    temporary.replace(args.allc_table)
    summary.update({
        "allc_built": built,
        "allc_reused_local": reused_local,
        "allc_reused_existing": reused_existing,
        "allc_table": str(args.allc_table),
    })
    (args.allc_table.parent / "prepare_summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
