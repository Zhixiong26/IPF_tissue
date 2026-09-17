# IPF tissue MethylVI workflows

This directory contains three coordinated MethylVI feature routes over the same
6,264 MethSCAn-filtered cells. The completed baseline experiment contains eight
ALLCools/VMR models; the VMR+pooled-DMR extension adds six models.

```text
Bismark coverage / ALLC
├── allcools/       ALLCools 5-kb feature selection → MethylVI
├── MethSCAn output VMRs
│   └── vmr/        selected MethSCAn VMR BED → MethylVI
└── pooled cell-type DMRs + selected VMRs
    └── vmr_dmr/    all unique pooled hypo-DMRs + VMR → MethylVI

shared/             training, supervised-UMAP, and methylation-QC code
```

Production one-click submission:

```bash
bash Scripts/Methylvi/submit_methylvi_8models.sh
```

The dependency graph builds shared inputs first, trains eight models in
parallel when their inputs are ready, and runs a final output summary. The
experiment root must not already exist.

The completed baseline result is `Results/MethylVI_clean6264_10k30k` (8/8
models complete). The VMR+DMR extension is submitted separately:

```bash
bash Scripts/Methylvi/vmr_dmr/submit_vmr_dmr_pipeline.sh
```

It writes only to `Results/MethylVI_vmr_plus_all_unique_pooled_DMR`, reuses the
completed baseline VMR inputs, and serializes its six cu03 join/train pairs.

Individual route runners remain available for verification and diagnosis:

```bash
bash Scripts/Methylvi/allcools/run.sh <stage>
bash Scripts/Methylvi/vmr/run.sh <stage>
```

All routes use `03_filtered/column_header.txt` as the sole cell-selection
authority. The selected-ALLC manifest maps those cell IDs to original indexed
ALLCs. The VMR route also requires a completed
MethSCAn `run_summary.json` before preparing counts. The current completed
upstream run is `Results/Methscan/CYL_ZCP_full_20260826_final`.

The eight models are ALLCools 5-kb × {10k,30k} and MethSCAn VMR
{0.01,0.02,0.05} × {10k,30k}. MethylVI uses integer CGN `mc/cov`, trains with
`sample_id` as batch key, and exports ordinary UMAP/Leiden, supervised UMAP,
sequencing-depth, overall-mCG, and mean-mCG diagnostics.

The VMR+DMR route uses all pooled hypo-DMRs passing raw `p < 0.01` and
`abs(methdiff) >= 0.25`, without Top200/Top-N truncation. Exact intervals are
deduplicated within hypo cell type, overlapping DMRs are merged, and only DMRs
with zero overlap to a model's selected VMRs are appended. Its preparation
produced 1,104,139 exact-unique and 321,910 merged DMR intervals.

Read the corresponding workflow README before running. Override the run,
variance branch, or BED only with a completed MethSCAn result.
