# IPF 30wcov MethylVI

`30wcov` contains 6,554 plain-text Bismark coverage files, not ALLCools MCDS files. The builder therefore aggregates columns 5 (`mc`) and 5+6 (`cov`) into 5-kb bins directly and writes the `mCG.layers['mc']` / `mCG.layers['cov']` input expected by MethylVI.

The supplied `Supplementary/hg38.canonical.chrom.sizes` is used by default to retain only chr1–22, chrX and chrY bins within their GRCh38 bounds.

Use a Python >=3.9 environment containing `scvi-tools` with `scvi.external.METHYLVI`, `anndata`, `mudata`, `scanpy`, `torch`, `igraph`, and `leidenalg`.

The reference project's ENCODE GRCh38 blacklist is included at `Supplementary/ENCFF356LFX_GRCh38_blacklist.bed.gz` and is enabled by default. The builder excludes a 5-kb bin when blacklist overlap exceeds 20% (change with `IPF_BLACKLIST_FRACTION`).

Run a small end-to-end validation first:

```bash
cd /home/lijia/luozhixiong/IPF_tissue/Scripts
source 00_methylvi_config.sh
python3 01_build_cov_methylvi_input.py --max-cells 100 --top-bins 5000 --min-cells 10 --threads 8 \
  --output "$IPF_MVI_ROOT/smoke/input.h5mu" --work-dir "$IPF_MVI_ROOT/smoke/checkpoints"
python3 02_train_methylvi.py --input "$IPF_MVI_ROOT/smoke/input.h5mu" --output "$IPF_MVI_ROOT/smoke/results" --epochs 2 --accelerator cpu
```

Then run `bash 03_run_methylvi.sh all`. The default feature selection retains the 50,000 most variable bins observed at >=20 reads in >=50 cells. This is a direct-data substitute for the reference project’s ALLCools hypo-bin selection, not an identical feature-selection algorithm.

The supplied annotation has only `manual_celltype`; filenames provide two prefixes (`CYL`, `ZCP`), which the scripts use as `cohort` and MethylVI batch covariate. Add donor/disease metadata before interpreting disease effects, and set `IPF_BATCH_KEY` to the appropriate technical batch column.
