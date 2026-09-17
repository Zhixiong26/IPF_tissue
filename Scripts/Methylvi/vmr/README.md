# MethSCAn VMRs → MethylVI

This workflow trains MethylVI on integer CGN `mc/cov` counts aggregated over a
VMR BED produced by a completed MethSCAn run. It is independent from the
ALLCools 5-kb workflow and writes to a separate result root.

## Choose the MethSCAn VMR set

MethSCAn produces one VMR BED per variance threshold. The defaults target the
current formal run `CYL_ZCP_full_20260826_final` and its `0.01` branch. To use another completed run or branch, override the run/threshold or
the BED directly; for example:

```bash
export IPF_METHSCAN_VMR_SOURCE="/home/lijia/luozhixiong/IPF_tissue/Results/Methscan/<run>/04_scan/var_0.01/VMRs.bed"
```

`run.sh verify` refuses to proceed when the selected BED is missing. The VMR
builder uses `03_filtered/column_header.txt` as the sole cell list and uses
the selected-ALLC manifest only to resolve paths, so MethylVI uses the same 6,264 post-filter
cells. VMRs are restricted to canonical chromosomes, blacklist-filtered,
checked for overlap, and assigned the highest covered-cell cutoff that still
permits the 30k target. The resulting cell cutoff and effective percentage are
recorded in `build_summary.json`.

After coverage eligibility, VMRs are ranked by `peak_var` descending, then
`n_obs_cells`, `n_cpg`, and stable VMR ID. Each threshold builds exact top 30k
once and derives a nested top 10k without rescanning ALLCs.

## Stages

```bash
bash Scripts/Methylvi/vmr/run.sh verify
bash Scripts/Methylvi/vmr/run.sh prepare
bash Scripts/Methylvi/vmr/run.sh build
bash Scripts/Methylvi/vmr/run.sh train
bash Scripts/Methylvi/vmr/run.sh plots
bash Scripts/Methylvi/vmr/run.sh all
```

Run all three thresholds at both feature targets through the project-level DAG:

```bash
bash Scripts/Methylvi/submit_methylvi_8models.sh
```

The default output root is `Results/MethylVI_30wcov_vmrs_blacklist_f0p2`.
Paths and MethylVI parameters are in `00_vmr_methylvi_config.sh`. Submit the
full cluster job with:

```bash
sbatch Scripts/Methylvi/vmr/slurm/run_methylvi_vmrs.sbatch
```

## Formal completed experiment

`Results/MethylVI_clean6264_10k30k` contains six completed VMR models:
`{0.01,0.02,0.05} x {Top10k,Top30k}`. Every route has `model.COMPLETE`, a
20-dimensional latent embedding, and UMAPs colored by cell type, sample,
condition, and MethylVI Leiden cluster. The same selected VMR inputs are reused
by `../vmr_dmr/`, which appends only non-overlapping pooled DMR features and
writes to a separate result root.
