# Shared MethylVI steps

`04_train_methylvi.py`, `05_plot_supervised_umap.py`, and
`07_plot_methylation_qc.py` are shared by the ALLCools, MethSCAn-VMR, and
VMR+pooled-DMR workflows. They consume common `IPF_*` MethylVI environment
variables, which each workflow runner sets for its own result route.

Training uses the scvi-tools default random cell split with seed 0: for the
6,264-cell production inputs this is 5,638 training cells, 626 validation
cells, and no independent test set. The validation ELBO drives early stopping;
latent embeddings and UMAPs are subsequently generated for all 6,264 cells.
