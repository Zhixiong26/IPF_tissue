#!/usr/bin/env python3
"""Calculate mCG, mCH, mCCC and covered CpG sites from per-cell ALLC files."""

import argparse
import concurrent.futures as cf
import csv
from pathlib import Path

from qc_common import atomic_tsv, open_text, value_or_na


FIELDS = ["sample_id", "barcode", "cell_id", "mc_CG", "cov_CG", "mCG",
          "mc_CH", "cov_CH", "mCH", "mc_CCC", "cov_CCC", "mCCC",
          "n_CpG_covered", "allc_context_status", "allc_path"]


def inputs(path):
    with Path(path).open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or not {"sample_id", "cell_id", "allc_path"}.issubset(rows[0]):
        raise ValueError("ALLC manifest needs sample_id, cell_id, allc_path: %s" % path)
    return rows


def ratio(mc, cov):
    return mc / cov if cov else None


def process_allc(item):
    sample_id, cell_id, path = item["sample_id"], item["cell_id"], Path(item["allc_path"])
    barcode = item.get("barcode") or (cell_id[len(sample_id) + 1:] if cell_id.startswith(sample_id + "_") else cell_id)
    mc_cg = cov_cg = mc_ch = cov_ch = mc_ccc = cov_ccc = n_cpg = 0
    contexts = set()
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                raise ValueError("%s:%s has fewer than 6 columns" % (path, line_number))
            context = fields[3].upper()
            mc, cov = int(fields[4]), int(fields[5])
            if mc < 0 or cov < mc:
                raise ValueError("Invalid mc/cov at %s:%s" % (path, line_number))
            contexts.add(context)
            if context.startswith("CG"):
                mc_cg += mc
                cov_cg += cov
                if cov > 0:
                    n_cpg += 1
            elif len(context) >= 2 and context[0] == "C" and context[1] != "G":
                mc_ch += mc
                cov_ch += cov
                if context == "CCC":
                    mc_ccc += mc
                    cov_ccc += cov
    has_ch = any(len(c) >= 2 and c.startswith("C") and c[1] != "G" for c in contexts)
    has_ccc = "CCC" in contexts
    status = "full_CG_CH_CCC" if has_ch and has_ccc else ("CCC_absent" if has_ch else "CG_only_CH_CCC_unavailable")
    return {"sample_id": sample_id, "barcode": barcode, "cell_id": cell_id,
            "mc_CG": mc_cg, "cov_CG": cov_cg, "mCG": value_or_na(ratio(mc_cg, cov_cg)),
            "mc_CH": value_or_na(mc_ch if has_ch else None),
            "cov_CH": value_or_na(cov_ch if has_ch else None),
            "mCH": value_or_na(ratio(mc_ch, cov_ch) if has_ch else None),
            "mc_CCC": value_or_na(mc_ccc if has_ccc else None),
            "cov_CCC": value_or_na(cov_ccc if has_ccc else None),
            "mCCC": value_or_na(ratio(mc_ccc, cov_ccc) if has_ccc else None),
            "n_CpG_covered": n_cpg, "allc_context_status": status,
            "allc_path": str(path.resolve())}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=16,
                        help="Parallel per-cell ALLC files; keep modest on shared storage")
    parser.add_argument("--max-files", type=int, default=0, help="Testing only; 0 scans all files")
    args = parser.parse_args()
    if args.workers < 1 or args.max_files < 0:
        raise ValueError("--workers must be >0 and --max-files must be >=0")
    tasks = inputs(args.manifest)
    if args.max_files:
        tasks = tasks[:args.max_files]
    effective_workers = min(args.workers, len(tasks))
    rows = []
    if effective_workers == 1:
        iterator = map(process_allc, tasks)
        executor = None
    else:
        executor = cf.ProcessPoolExecutor(max_workers=effective_workers)
        iterator = executor.map(process_allc, tasks)
    try:
        for number, row in enumerate(iterator, start=1):
            rows.append(row)
            if number % 500 == 0 or number == len(tasks):
                print("ALLC complete: %s/%s" % (number, len(tasks)), flush=True)
    finally:
        if executor is not None:
            executor.shutdown()
    if not rows:
        raise RuntimeError("No ALLC files were processed")
    atomic_tsv(args.output, FIELDS, rows)


if __name__ == "__main__":
    main()
