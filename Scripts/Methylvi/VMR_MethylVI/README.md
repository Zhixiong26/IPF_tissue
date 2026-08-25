# IPF VMR-MethylVI

This pipeline trains MethylVI on integer mCG counts aggregated over the IPF VMR set. It is independent of the 5-kb analysis and does not overwrite its outputs.

## Inputs and filtering

- Source regions: `Data/cov/VMR_1%.txt` (37,798 intervals).
- Only chromosomes in `Supplementary/hg38.canonical.chrom.sizes` are retained.
- VMRs with at least 20% overlap with `ENCFF356LFX_GRCh38_blacklist.bed.gz` are removed.
- The prepared BED contains 36,025 non-overlapping VMRs with unique IDs.
- During count assembly, a VMR is retained when `cov > 0` in more than 200 of 6,554 cells.
- `CYL`/`ZCP` (`cohort`) is used as the MethylVI batch.

## Data model

Each ALLC is scanned once. For every retained VMR, the pipeline sums integer methylated counts (`mc`) and total coverage (`cov`) in CGN context. These matrices are stored in `mCG.layers['mc']` and `mCG.layers['cov']` of a MuData file and passed directly to MethylVI.

The old `Data/cov/30_40_VMRs.sh` Hamming-distance workflow is not used because binary CpG calls and pairwise Hamming distances are not MethylVI inputs.

## Commands

```bash
cd /home/lijia/luozhixiong/IPF_tissue
bash VMR_MethylVI/Scripts/04_run_vmr_methylvi.sh verify
bash VMR_MethylVI/Scripts/04_run_vmr_methylvi.sh prepare
bash VMR_MethylVI/Scripts/04_run_vmr_methylvi.sh build
bash VMR_MethylVI/Scripts/04_run_vmr_methylvi.sh train
bash VMR_MethylVI/Scripts/04_run_vmr_methylvi.sh plots
bash VMR_MethylVI/Scripts/04_run_vmr_methylvi.sh supervised
```

Full Slurm run (not submitted automatically):

```bash
sbatch VMR_MethylVI/run_methylvi_vmrs.sbatch
```

Default output root:

```text
Results/MethylVI_30wcov_vmrs_blacklist_f0p2
```
