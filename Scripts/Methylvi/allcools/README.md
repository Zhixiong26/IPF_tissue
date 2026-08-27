# ALLCools → MethylVI

This workflow starts with Bismark coverage/ALLC files, uses ALLCools 5-kb
hypo-score regions for feature selection and clustering, then reconstructs
integer CGN `mc/cov` counts for MethylVI. The ALLCools score itself is never
used as a MethylVI count matrix.

## Stages

| Command | Result |
|---|---|
| `verify` | Validate input paths, environments, and cell inventory. |
| `prepare` | Reuse or create ALLC files and generate the 5-kb MCDS. |
| `cluster` | Run ALLCools filtering, LSI, and clustering. |
| `build` | Aggregate retained bins to integer MethylVI `mc/cov` layers. |
| `train` | Train MethylVI and save the latent embedding. |
| `plots-before` / `plots-after` / `supervised` | Export the respective UMAPs. |

Run from the repository root:

```bash
bash Scripts/Methylvi/allcools/run.sh verify
bash Scripts/Methylvi/allcools/run.sh all
```

The default output root is
`Results/MethylVI_30wcov_allcools_blacklist_f0p2`. Configure paths, feature
selection, environments, and training parameters in `00_methylvi_config.sh`
or by exporting the documented `IPF_*` variables before invoking `run.sh`.

The supplied Slurm jobs are in `slurm/`; submit them from the repository root,
for example `sbatch Scripts/Methylvi/allcools/slurm/run_methylvi_allcools.sbatch`.

