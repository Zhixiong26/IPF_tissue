# MethSCAn VMRs → MethylVI

This workflow trains MethylVI on integer CGN `mc/cov` counts aggregated over a
VMR BED produced by a completed MethSCAn run. It is independent from the
ALLCools 5-kb workflow and writes to a separate result root.

## Choose the MethSCAn VMR set

MethSCAn produces one VMR BED per variance threshold. The defaults target the
current formal run `CYL_ZCP_full_20260826_final` and its `0.01` branch. That
run is not yet scan-complete, so `verify` deliberately stops until its VMR BED
exists. To use another completed run or branch, override the run/threshold or
the BED directly; for example:

```bash
export IPF_METHSCAN_VMR_SOURCE="/home/lijia/luozhixiong/IPF_tissue/Results/Methscan/<run>/04_scan/var_0.01/VMRs.bed"
```

`run.sh verify` refuses to proceed when the selected BED is missing. The VMR
builder reads the paired MethSCAn selected-ALLC manifest, preserving the exact
upstream cell set and original indexed ALLC paths. VMRs are restricted to
canonical chromosomes, blacklist-filtered, checked for overlap, and retained
only when covered in more than the configured number of cells.

## Stages

```bash
bash Scripts/Methylvi/vmr/run.sh verify
bash Scripts/Methylvi/vmr/run.sh prepare
bash Scripts/Methylvi/vmr/run.sh build
bash Scripts/Methylvi/vmr/run.sh train
bash Scripts/Methylvi/vmr/run.sh plots
bash Scripts/Methylvi/vmr/run.sh all
```

The default output root is `Results/MethylVI_30wcov_vmrs_blacklist_f0p2`.
Paths and MethylVI parameters are in `00_vmr_methylvi_config.sh`. Submit the
full cluster job with:

```bash
sbatch Scripts/Methylvi/vmr/slurm/run_methylvi_vmrs.sbatch
```
