# IPF tissue MethylVI workflows

This directory contains two independent, reproducible MethylVI analysis lines.

```text
Bismark coverage / ALLC
├── allcools/       ALLCools 5-kb feature selection → MethylVI
└── MethSCAn output VMRs
    └── vmr/        selected MethSCAn VMR BED → MethylVI

shared/             training and supervised-UMAP code used by both lines
```

Use only one workflow runner at a time:

```bash
bash Scripts/Methylvi/allcools/run.sh <stage>
bash Scripts/Methylvi/vmr/run.sh <stage>
```

Read the corresponding workflow README before running. The VMR route defaults
to the latest MethSCAn run's selected-ALLC manifest, but it remains blocked
until that run has produced the requested `VMRs.bed` branch. Override the run,
variance branch, or BED only with a completed MethSCAn result.
