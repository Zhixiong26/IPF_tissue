#!/usr/bin/env python3
"""Outer-join FASTQ and BAM per-cell counts and calculate mapping efficiency."""

import argparse
from pathlib import Path

from qc_common import atomic_tsv, cell_key, read_tsv, required_columns, value_or_na


FIELDS = ["sample_id", "barcode", "cell_id", "input_pairs", "input_reads",
          "mapped_pairs", "mapped_reads", "mapping_rate", "final_reads",
          "final_reads_definition", "count_status"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fastq", type=Path, action="append", required=True,
                        help="Per-sample step-01 table; repeat to combine samples")
    parser.add_argument("--bam", type=Path, action="append", required=True,
                        help="Per-sample step-02 table; repeat to combine samples")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    fastq_rows = [row for path in args.fastq for row in read_tsv(path)]
    bam_rows = [row for path in args.bam for row in read_tsv(path)]
    required_columns(fastq_rows, ("sample_id", "barcode", "input_pairs", "input_reads"), args.fastq)
    required_columns(bam_rows, ("sample_id", "barcode", "mapped_pairs", "mapped_reads"), args.bam)
    fastq = {cell_key(row): row for row in fastq_rows}
    bam = {cell_key(row): row for row in bam_rows}
    if len(fastq) != len(fastq_rows) or len(bam) != len(bam_rows):
        raise ValueError("Duplicate sample_id/barcode key in input table")
    output = []
    for sample_id, barcode in sorted(set(fastq) | set(bam)):
        fq, alignment = fastq.get((sample_id, barcode)), bam.get((sample_id, barcode))
        input_pairs = int(fq["input_pairs"]) if fq else None
        input_reads = int(fq["input_reads"]) if fq else None
        mapped_pairs = int(alignment["mapped_pairs"]) if alignment else 0
        mapped_reads = int(alignment["mapped_reads"]) if alignment else 0
        rate = mapped_pairs / input_pairs if input_pairs else None
        status = "complete" if fq and alignment else ("missing_bam" if fq else "missing_fastq")
        output.append({"sample_id": sample_id, "barcode": barcode,
                       "cell_id": "%s_%s" % (sample_id, barcode),
                       "input_pairs": value_or_na(input_pairs), "input_reads": value_or_na(input_reads),
                       "mapped_pairs": mapped_pairs, "mapped_reads": mapped_reads,
                       "mapping_rate": value_or_na(rate), "final_reads": mapped_reads,
                       "final_reads_definition": "primary_mapped_records_R1_plus_R2",
                       "count_status": status})
    atomic_tsv(args.output, FIELDS, output)


if __name__ == "__main__":
    main()
