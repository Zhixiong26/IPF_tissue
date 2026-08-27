# Shared MethylVI steps

`04_train_methylvi.py` and `05_plot_supervised_umap.py` are shared by the
ALLCools and MethSCAn-VMR workflows. They consume the common `IPF_*` MethylVI
environment variables, which each workflow runner sets from its own config.

