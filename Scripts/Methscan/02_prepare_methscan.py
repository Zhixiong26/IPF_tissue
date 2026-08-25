#!/usr/bin/env python3
"""Run `methscan prepare --input-format allc` from the canonical manifest."""

import argparse
import csv
import json
import subprocess
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--methscan", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--chunksize", type=int, default=10000000)
    parser.add_argument("--round-sites", action="store_true")
    args = parser.parse_args()
    if not args.methscan.is_file():
        raise FileNotFoundError(args.methscan)
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise FileExistsError("Prepare output is not empty: %s" % args.output_dir)
    with args.manifest.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or "link_path" not in rows[0]:
        raise ValueError("Manifest has no input rows or link_path column")
    paths = [row["link_path"] for row in rows]
    missing = [path for path in paths if not Path(path).is_file()]
    if missing:
        raise FileNotFoundError("Missing manifest inputs, first: %s" % missing[0])
    command = [
        str(args.methscan), "prepare", "--input-format", "allc",
        "--chunksize", str(args.chunksize),
    ]
    if args.round_sites:
        command.append("--round-sites")
    command.extend(paths)
    command.append(str(args.output_dir))
    args.output_dir.parent.mkdir(parents=True, exist_ok=True)
    provenance = {
        "command_prefix": command[:7],
        "input_count": len(paths),
        "manifest": str(args.manifest.resolve()),
        "output_dir": str(args.output_dir.resolve()),
        "round_sites": args.round_sites,
    }
    with (args.output_dir.parent / "prepare_command.json").open("w") as handle:
        json.dump(provenance, handle, indent=2, sort_keys=True)
        handle.write("\n")
    subprocess.run(command, check=True)


if __name__ == "__main__":
    main()
