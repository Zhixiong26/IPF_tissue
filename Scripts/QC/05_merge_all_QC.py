#!/usr/bin/env python3
"""Merge sequencing/mapping and ALLC metrics, retaining every QC flag."""

import argparse
from pathlib import Path

from qc_common import atomic_tsv, cell_key, parse_float, parse_int, read_tsv, required_columns


def flag(value, predicate):
    return "NA" if value is None else str(bool(predicate(value))).upper()


def bool_value(value):
    return None if value == "NA" else value == "TRUE"


def conjunction(values):
    values = [bool_value(value) for value in values]
    return "NA" if any(value is None for value in values) else str(all(values)).upper()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts", type=Path, required=True)
    parser.add_argument("--allc", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--min-final-reads", type=int, default=500000)
    parser.add_argument("--max-final-reads", type=int, default=10000000)
    parser.add_argument("--min-mapping-rate", type=float, default=0.50)
    parser.add_argument("--min-mcg", type=float, default=0.50)
    parser.add_argument("--max-mch", type=float, default=0.20)
    parser.add_argument("--max-mccc", type=float, default=0.05)
    parser.add_argument("--min-cpg", type=int, default=300000)
    parser.add_argument("--max-cpg", type=int, default=1200000)
    args = parser.parse_args()
    counts, allc = read_tsv(args.counts), read_tsv(args.allc)
    required_columns(counts, ("sample_id", "barcode", "mapping_rate", "final_reads"), args.counts)
    required_columns(allc, ("sample_id", "barcode", "mCG", "mCH", "mCCC", "n_CpG_covered"), args.allc)
    left, right = {cell_key(r): r for r in counts}, {cell_key(r): r for r in allc}
    if len(left) != len(counts) or len(right) != len(allc):
        raise ValueError("Duplicate sample_id/barcode key")
    fields = ["sample_id", "barcode", "cell_id", "input_pairs", "input_reads", "mapped_pairs",
              "mapped_reads", "mapping_rate", "final_reads", "final_reads_definition",
              "mCG", "mCH", "mCCC", "n_CpG_covered", "allc_context_status",
              "pass_reads", "pass_mapping", "pass_mCG", "pass_mCH", "pass_mCCC",
              "pass_primary_qc", "pass_CpG", "pass_final_qc", "missing_data", "fail_reasons"]
    output = []
    for key in sorted(set(left) | set(right)):
        count, methyl = left.get(key, {}), right.get(key, {})
        final_reads, mapping = parse_int(count.get("final_reads")), parse_float(count.get("mapping_rate"))
        mcg, mch, mccc = parse_float(methyl.get("mCG")), parse_float(methyl.get("mCH")), parse_float(methyl.get("mCCC"))
        cpg = parse_int(methyl.get("n_CpG_covered"))
        flags = {
            "pass_reads": flag(final_reads, lambda x: args.min_final_reads < x < args.max_final_reads),
            "pass_mapping": flag(mapping, lambda x: x > args.min_mapping_rate),
            "pass_mCG": flag(mcg, lambda x: x > args.min_mcg),
            "pass_mCH": flag(mch, lambda x: x < args.max_mch),
            "pass_mCCC": flag(mccc, lambda x: x < args.max_mccc),
            "pass_CpG": flag(cpg, lambda x: args.min_cpg <= x <= args.max_cpg),
        }
        flags["pass_primary_qc"] = conjunction([flags[x] for x in ("pass_reads", "pass_mapping", "pass_mCG", "pass_mCH", "pass_mCCC")])
        flags["pass_final_qc"] = conjunction([flags["pass_primary_qc"], flags["pass_CpG"]])
        missing = [name for name, value in (("FASTQ", count.get("input_pairs")), ("mapping_rate", count.get("mapping_rate")),
                                             ("ALLC", methyl.get("mCG")), ("mCH", methyl.get("mCH")),
                                             ("mCCC", methyl.get("mCCC"))) if value is None or value in ("", "NA")]
        failures = [name.replace("pass_", "") for name, value in flags.items() if name not in ("pass_primary_qc", "pass_final_qc") and value == "FALSE"]
        row = {name: "NA" for name in fields}
        for source in (count, methyl):
            for name in fields:
                if name in source:
                    row[name] = source[name]
        row.update(flags)
        row.update({"sample_id": key[0], "barcode": key[1], "cell_id": "%s_%s" % key,
                    "missing_data": ";".join(missing), "fail_reasons": ";".join(failures)})
        output.append(row)
    atomic_tsv(args.output, fields, output)


if __name__ == "__main__":
    main()
