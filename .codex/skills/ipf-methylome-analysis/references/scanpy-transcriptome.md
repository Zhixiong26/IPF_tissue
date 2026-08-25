# CYL/ZCP transcriptome-E Scanpy workflow

Read this reference whenever the task concerns RNA/transcriptome-E QC, Scrublet, Scanpy, Harmony, Leiden clustering, marker review, cell-type annotation, result export, or Slurm submission.

## Scope and source of truth

- Project root: `/home/lijia/luozhixiong/IPF_tissue`
- Included cohorts: CYL and ZCP only. LC is not part of the current transcriptome analysis.
- Inputs:
  - `Data/Matrix/25100718_CYL_E/filtered_feature_bc_matrix.zip`
  - `Data/Matrix/25100718_ZCP_E/filtered_feature_bc_matrix.zip`
- Canonical and preferred execution entry point: `Scripts/Scanpy/Notebooks/E_CYL_ZCP_scanpy.ipynb`
- Secondary batch conversion: `Scripts/Scanpy/Scripts/01_run_e_scanpy.py`
- Secondary Slurm wrapper: `Scripts/Scanpy/Scripts/run_e_scanpy.sbatch`
- Operational guide: `Scripts/Scanpy/README.md`
- Complete parameter/evidence record: `Scripts/Scanpy/Report.md`

Treat `E_CYL_ZCP_scanpy.ipynb` and the stage Report as the maintained source of truth. Default Scanpy execution to a clean top-to-bottom replay of that notebook in a new result root. Do not use `Scripts/Scanpy/Scripts/` merely because a Python/Slurm entry point exists: it is secondary, and the current conversion has already produced a divergent 17-cluster result. Use it only at explicit user request, for validation/export helpers, or after equivalence to the notebook has been demonstrated for retained cell IDs, Leiden assignments, markers, and cell types in the intended execution environment. Do not copy parameters from the older `S12-2N.ipynb` or the pre-validation draft.

## Execution priority

1. Run the canonical notebook with the `ipf-allcools` kernel, from the first cell, into a new named `Results/Scanpy/` replay root. Never overwrite the source notebook.
2. Verify input/QC/singlet counts, observed Leiden set, ranked markers, annotation guard, figures, and expected notebook errors such as an intentionally closed save gate.
3. After the user confirms the reviewed notebook result, export its H5AD, per-cell cell-type table, plots, parameters, and replay summary from the same notebook-derived state.
4. Treat `Scripts/Scanpy/Scripts/` as secondary. Before promoting a script result, compare its retained cell IDs and per-cell Leiden/cell-type assignments with the canonical notebook result; matching aggregate counts and parameters alone are insufficient.

For annotation evidence hierarchy, the reviewed decisions for clusters 0–17, and the required correction procedure, also read `scanpy-cluster-validation.md`.

## Confirmed workflow contract

1. Read CYL and ZCP separately, prefix cell IDs with cohort, and retain integer counts in `layers['counts']`.
2. Calculate mitochondrial and ribosomal QC metrics per cohort.
3. Assess and filter each cohort independently using the recorded QC thresholds.
4. Run Scrublet independently on each filtered cohort using raw counts and the recorded dynamic expected-doublet-rate formula.
5. Merge retained singlets with an outer gene join.
6. Normalize to 10,000 counts per cell, apply `log1p`, and preserve the full normalized gene matrix in `raw`.
7. Select cohort-aware HVGs with the recorded Seurat flavor; do not use `seurat_v3` in the current environment.
8. Scale HVGs without covariate regression, calculate PCA, and run Harmony with `cohort` as the batch key.
9. Build matched neighbor graphs and UMAPs from original `X_pca` and corrected `X_pca_harmony`; calculate Leiden only on the Harmony graph.
10. Rank cluster markers with Wilcoxon on `raw`, generate the canonical-marker dotplot, and apply the reviewed annotation only when the observed cluster ID set exactly matches the reviewed set.
11. Export machine-readable QC, doublet, HVG, PCA, cluster, marker, annotation, and per-cell tables alongside figures and H5AD.

The exact defaults—including QC cutoffs, PCA/Harmony/KNN/Leiden settings, marker panels, cluster annotation evidence, plot settings, output schema, and Slurm resources—are recorded in `Scripts/Scanpy/Report.md`. Read that file before changing parameters, interpreting results, or preparing a submission.

## Annotation invariant

The currently reviewed mapping covers Leiden clusters 0–17 from the saved notebook result. Cluster 15 has combined lymphatic/pan-endothelial evidence; cluster 14 remains a cycling label, and user-specified clusters 16/17 share the primary label `NA` while retaining their distinct low-confidence evidence notes. Cluster numbers are not biological identifiers and may change after filtering, dependency, or parameter changes.

- Never silently apply the 0–17 mapping to a different observed cluster set.
- The production script must retain the cluster-set guard.
- When IDs differ, retain diagnostic marker tables, plots, and an unreviewed H5AD; label every cell `Unassigned` because even retained numeric IDs may represent different biology, and require marker review before updating the mapping.
- `--allow-unreviewed-clusters` changes process exit behavior only. It does not make the annotation valid.
- Record annotation evidence and confidence in both the stage Report and exported annotation table.
- Generate the focused epithelial panel for clusters 1/4/6/10/12/14/16 and the lineage/QC audit for rare clusters 15/17 before promoting a provisional subtype label.
- Treat the prior broad CYL/ZCP annotation as lineage-level auxiliary evidence, not ground truth. Do not replace current T-cell or lymphatic-endothelial calls with prior `NA`; a cell-level comparison requires barcode-resolved prior labels.

## Output and plotting invariant

- Immutable/local inputs belong under `Data/`; derived analysis outputs belong under a new named root in `Results/Scanpy/`.
- For the canonical CYL/ZCP notebook, every figure displayed by a plotting cell must also be saved as a 300-dpi PNG under `Results/Scanpy/E_CYL_ZCP_notebook/figures/`, and the final cell must write `figure_manifest.json` with the exact saved paths and count.
- A figure-only notebook replay must not overwrite an existing confirmed H5AD, per-cell annotation table, or parameter JSON. Preserve them only after exact per-cell annotation equivalence is verified; otherwise stop for review. Explicitly enable data overwrite only when the user intends to replace formal results.
- Scheduler stdout/stderr belong in `Scripts/Scanpy/logs/`.
- Refuse non-empty formal result roots. Use a distinct run name rather than overwriting a prior run.
- Save a sample-colored before/after-Harmony UMAP comparison using matched parameters.
- Cell-type, sample, and group UMAPs are separate Harmony-coordinate figures with legends in the right margin. Do not put long cell-type labels on the embedding.
- If the cluster-ID guard fires, warn the user that the original run intentionally lacks a final cell-type UMAP. Once every changed cluster has marker-supported review and a guarded derivative H5AD, a post-run cell-type UMAP may be generated from the preserved Harmony coordinates; label it as a reviewed derivative and save the warning/status alongside the figure.
- Current sample values are CYL and ZCP. No verified biological group is available, so skip the group figure with an explicit status record rather than inventing a group from sample identity.
- Preserve tabular evidence for every plotted result so the analysis does not depend on visual inspection alone.
- Merge reviewed clusters by `cell_type` before marker-dotplot aggregation. Use exactly the same cell-type labels and dendrogram-synchronized order for the y-axis and top marker groups; show three to four available markers per type, red mean-expression color, fraction-sized dots, and the right-side dendrogram.

The final H5AD uses scaled HVGs in `X`, the full log-normalized gene matrix in `raw`, retained raw counts for the HVG work object in `layers['counts']`, Harmony/PCA/UMAP representations in `obsm`, and QC/doublet/cohort/cluster/cell-type metadata in `obs`.

## Environment and execution

- Verified Python: `/home/lijia/jiangyuanpei/miniforge3/envs/allcools/bin/python`
- Confirmed baseline: Python 3.9, Scanpy 1.9.3, harmonypy 0.0.10, Scrublet 0.2.3.
- Do not switch to the Python 3.12 `methylvi` environment for this workflow.
- A Jupyter kernel failure caused by Anaconda channel Terms of Service is an environment-launch issue, not a reason to mutate the shared environment. Prefer the working allcools interpreter and obtain approval before environment changes.
- Before a notebook replay, verify the `ipf-allcools` kernel, both matrix ZIPs, a new result root, and available storage. Preserve the executed notebook and a machine-readable replay summary.
- Before an explicitly requested Slurm/script run, run Python syntax/import checks, `bash -n` on the wrapper, confirm both matrix ZIPs are readable, ensure the result root is new/empty, and check available storage.
- Do not submit a production script merely because it exists. Notebook execution remains the default; use the script only under the secondary-entry conditions above.

## Compute-node path boundary

The fat compute node cannot directly access `/mnt/data04`. A project symlink whose target is on `/mnt/data04` is therefore not a valid compute-node input. Resolve and verify actual paths before submission.

For the current transcriptome workflow, use the copied ZIP matrices under the project `Data/Matrix/` paths above and write outputs under project `Results/Scanpy/`. Do not redirect the job to `/mnt/data04`.

## Completion evidence

A completed notebook execution or submitted job is not by itself a completed analysis. Before marking the Scanpy stage complete, verify:

- executed-notebook errors and save-gate status; for a script/Slurm run, also verify Slurm exit status, stderr, and `run_summary.json`;
- per-cohort input, QC, doublet, and retained-singlet counts;
- PCA/Harmony/UMAP outputs and cohort mixing diagnostics;
- observed Leiden cluster set;
- marker tables and annotation guard status;
- final H5AD readability and expected `obs`, `var`, `layers`, `raw`, `obsm`, and `uns` fields;
- expected figures and machine-readable tables are non-empty.

Record the command, job ID, result root, verified counts, annotation status, failures, and remaining decisions in `references/analysis-state.md` and `Scripts/Scanpy/Report.md`.
