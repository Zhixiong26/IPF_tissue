#!/usr/bin/env python3
"""Summarize per-cell QC flags by sample without deleting any cells."""

import argparse
from collections import defaultdict
from pathlib import Path

from qc_common import atomic_tsv, read_tsv, required_columns


FLAGS = ["pass_reads", "pass_mapping", "pass_mCG", "pass_mCH", "pass_mCCC",
         "pass_CpG", "pass_primary_qc", "pass_final_qc"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = read_tsv(args.input)
    required_columns(rows, ["sample_id"] + FLAGS, args.input)
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["sample_id"]].append(row)
    fields = ["sample_id", "total_cells", "fail_reads", "fail_mapping", "fail_mCG",
              "fail_mCH", "fail_mCCC", "fail_CpG", "missing_primary",
              "pass_primary", "pass_final"]
    output = []
    for sample_id, sample_rows in sorted(grouped.items()):
        output.append({"sample_id": sample_id, "total_cells": len(sample_rows),
                       "fail_reads": sum(r["pass_reads"] == "FALSE" for r in sample_rows),
                       "fail_mapping": sum(r["pass_mapping"] == "FALSE" for r in sample_rows),
                       "fail_mCG": sum(r["pass_mCG"] == "FALSE" for r in sample_rows),
                       "fail_mCH": sum(r["pass_mCH"] == "FALSE" for r in sample_rows),
                       "fail_mCCC": sum(r["pass_mCCC"] == "FALSE" for r in sample_rows),
                       "fail_CpG": sum(r["pass_CpG"] == "FALSE" for r in sample_rows),
                       "missing_primary": sum(r["pass_primary_qc"] == "NA" for r in sample_rows),
                       "pass_primary": sum(r["pass_primary_qc"] == "TRUE" for r in sample_rows),
                       "pass_final": sum(r["pass_final_qc"] == "TRUE" for r in sample_rows)})
    atomic_tsv(args.output, fields, output)


if __name__ == "__main__":
    main()
