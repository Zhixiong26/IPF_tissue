# ALLCools → MethylVI

This workflow starts from the 6,264 cells retained in the completed MethSCAn
filter header, resolves their original indexed ALLCs through the MethSCAn
input manifest, and uses ALLCools 5-kb
hypo-score regions for feature selection and clustering, then reconstructs
integer CGN `mc/cov` counts for MethylVI. The ALLCools score itself is never
used as a MethylVI count matrix.

## Stages

| Command | Result |
|---|---|
| `verify` | Validate input paths, environments, and cell inventory. |
| `prepare` | Reuse or create ALLC files and generate the 5-kb MCDS. |
| `cluster` | Rank blacklist-filtered 5-kb bins by current-cell hypo prevalence and select exact top 30k. |
| `build` | Aggregate retained bins to integer MethylVI `mc/cov` layers. |
| `train` | Train MethylVI and save the latent embedding. |
| `plots-before` / `plots-after` / `supervised` | Export the respective UMAPs. |

Run from the repository root:

```bash
bash Scripts/Methylvi/allcools/run.sh verify
bash Scripts/Methylvi/allcools/run.sh all
```

The production DAG builds 30k once and derives a strictly nested top-10k H5MU
using `selection_rank`; it does not rescan ALLCs. The effective
`MVI_HYPO_PERCENT` and hard feature-count check are written to
`feature_filter_summary.json`. Configure paths, feature
selection, environments, and training parameters in `00_methylvi_config.sh`
or by exporting the documented `IPF_*` variables before invoking `run.sh`.

The supplied Slurm jobs are in `slurm/`; submit them from the repository root,
for example `sbatch Scripts/Methylvi/allcools/slurm/run_methylvi_allcools.sbatch`.

The formal `Results/MethylVI_clean6264_10k30k` experiment has completed both
ALLCools models (Top10k and Top30k), including embeddings, ordinary and
supervised UMAPs, methylation QC, and `model.COMPLETE` markers.
