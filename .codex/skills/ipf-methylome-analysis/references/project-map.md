# Project map

## Inputs and current layout

- Project root: `/home/lijia/luozhixiong/IPF_tissue`
- GitHub intake requirement: for a new project, ask the user to create/provide the SSH repository URL and default branch, and confirm SSH readiness. Never request secrets or assume permission to commit or push.
- Required project directories: `Scripts/`, `Results/`, and `Supplementary/`.
- Initial matching stage directories under both `Scripts/` and `Results/`: `Environment/`, `FastQ/`, `QC/`, `Scanpy/`, `Methscan/`, and `Methylvi/`. Extend both parents together for new major stages.
- Each stage has a dedicated log directory at `Scripts/<stage>/logs/`; new jobs must not use a shared project-root log directory.
- Each `Scripts/<stage>/` contains a bilingual Chinese/English `README.md` (workflow contract and entry points) and `Report.md` (run/QC/result record). The layout helper creates only missing documents and never overwrites maintained content.
- Root-level `log/` and `logs/` are legacy directories. Preserve historical files; all new job and analysis logs belong in `Scripts/<stage>/logs/`.
- Linked raw data: `Data/Raw_fastq/E/` and `Data/Raw_fastq/Met/`
- Linked BAM data: `Data/Bam -> /mnt/data04/jiangyuanpei/xunyin_260727/data/bam`, containing the same six named E/Met batches. This link is visible on the login node but is not a valid fat-node input because `/mnt/data04` is inaccessible there.
- Current project selection: the user has directed the workflow to start from `Data/Bam` and temporarily skip FASTQ merging and processing. This does not retroactively validate upstream FASTQ QC or alignment quality.
- For future analyses, ask whether FASTQ is needed. Require a Raw FASTQ path only when the user chooses FASTQ; otherwise request the selected downstream data type, path, and known provenance. `Data/Raw_fastq/` is a project link location, not an implicit source selection.
- Current linked raw-data batches: `25100718_CYL`, `25100718_LC_S9_1N`, and `25100718_ZCP`, each with an `E` and a `Met` directory.
- The raw directory names are treated as provisional metadata only. Derive exact samples/read pairs from the FASTQ filenames before analysis.
- Raw-data validation entry point: `Scripts/FastQ/01_validate_raw_fastq.py`; outputs are written to `Results/FastQ/raw_fastq_validation/`.
- The server system Python is 3.6.8. Standalone intake/validation scripts must remain compatible with it unless they explicitly use a recorded Conda environment.
- Environment discovery entry point: `Scripts/Environment/01_discover_environments.sh`; inventory output belongs under `Results/Environment/`, while the approved shared/per-stage mapping belongs in `Supplementary/environments.tsv`.
- Per-cell QC entry points: `Scripts/QC/01_FASTQ_per_cell_reads.py` through `06_QC_summary.py`. They retain FASTQ, BAM, and ALLC layers plus individual flags and do not delete failed cells.
- Current environment strategy is hybrid: BAM and Scanpy/ALLCools share the read-only `allcools` environment; Methscan and MethylVI use separate verified environments; FASTQ/Bismark is deferred.

## CYL/ZCP transcriptome-E Scanpy workflow

- Current cohort scope: CYL and ZCP only; LC is excluded.
- Local compute-readable inputs:
  - `Data/Matrix/25100718_CYL_E/filtered_feature_bc_matrix.zip`
  - `Data/Matrix/25100718_ZCP_E/filtered_feature_bc_matrix.zip`
- Preferred/canonical entry point: `Scripts/Scanpy/Notebooks/E_CYL_ZCP_scanpy.ipynb`
- Secondary batch conversion: `Scripts/Scanpy/Scripts/01_run_e_scanpy.py`
- Secondary Slurm wrapper: `Scripts/Scanpy/Scripts/run_e_scanpy.sbatch`
- Complete parameter record: `Scripts/Scanpy/Report.md`
- Default result root: `Results/Scanpy/E_CYL_ZCP`; formal runs should use a new named/dated subdirectory.
- Environment: `/home/lijia/jiangyuanpei/miniforge3/envs/allcools/bin/python` with the recorded Python 3.9/Scanpy 1.9.3 baseline.
- Workflow: per-cohort QC and Scrublet, merge singlets, normalize/log1p, cohort-aware Seurat HVGs, scale/PCA, Harmony by cohort, Harmony neighbors, UMAP, Leiden, Wilcoxon markers, guarded manual annotation, and tabular/figure/H5AD export.
- The current notebook result has a reviewed mapping for clusters 0–17. Clusters 16/17 share the user-specified primary label `NA`; their differing low-confidence evidence remains documented. A fresh result with any other cluster set must not inherit this mapping automatically.
- A clean notebook replay reproduced the 18-cluster result, whereas the current Slurm script result diverged to 17 clusters. Default future Scanpy execution to the notebook. Do not promote a script result without explicit user direction or per-cell equivalence to the notebook in the intended execution environment.
- Read `references/scanpy-transcriptome.md` before Scanpy execution or interpretation.

## Existing coverage-to-MethylVI workflow

`Scripts/00_methylvi_config.sh` defines the project paths and two environment roots:

- ALLCools: `/home/lijia/jiangyuanpei/miniforge3/envs/allcools`
- MethylVI: `/home/lijia/luozhixiong/miniconda3/envs/methylvi`

The current entry point is `Scripts/03_run_methylvi.sh`:

| Stage | Command mode | Main output |
|---|---|---|
| Verify | `verify` | input/environment checks |
| ALLC preparation | `prepare` | ALLC input and MCDS |
| ALLCools clustering | `cluster` | 5-kb clustered H5AD |
| Build counts | `build` | MethylVI H5MU with integer mc/cov |
| Train | `train` | MethylVI latent embedding and model |
| Smoke test | `smoke` | isolated reduced test run |

`Scripts/README.md` documents the existing 6,554-cell coverage analysis and its invariants. It is the local source of truth for the current ALLCools/MethylVI implementation.

The root-level scripts and pre-existing `Results/MethylVI_*` directories predate the staged layout and remain valid historical work. New work uses the matching stage subdirectories.

## Boundaries to verify before new upstream work

- GitHub `origin` is `git@github.com:Zhixiong26/IPF_tissue.git`; the working branch observed on 2026-08-24 is `codex/initial-upload`. Authentication succeeds only when explicitly selecting `~/.ssh/id_ed25519_github`; no matching `IdentityFile` SSH configuration was observed. Push authorization remains unconfirmed.
- Genome build/reference and the Bismark genome directory are not yet established for the newly linked FASTQs.
- The exact library chemistry and the biological meaning of `E` versus `Met` are not yet established. For Met reads, barcode layout is established: R1 bases 1–17 are CB, R1 bases 18–29 are UMI, and the 8-bp header suffix is the sample index.
- For production per-cell FASTQ counts, run BAM step 02 first and pass its `sample_id/barcode` output to step 01 as the whitelist. This preserves mapped cells with zero/low additional input support while excluding error-derived raw barcode strings.
- Existing converted ALLCs contain only `CGN`; mCH and mCCC require full-context methylation calls and must remain missing for these CpG-only files.
- The maintained project task list is `Supplementary/TODO.md`.
- MethSCAn 1.1.0 is verified at `/home/lijia/jiangyuanpei/miniforge3/envs/MethSCAn/bin/methscan`. The maintained stage starts from indexed ALLCools ALLCs and performs native ALLC prepare, filter, smooth, VMR scan, matrix construction, and downstream Scanpy embedding/clustering through `Scripts/Methscan/run_methscan.sbatch`.
- Compute-node storage boundary: do not pass `/mnt/data04` or a symlink resolving there to fat-node jobs. Verify local staged inputs and free space before submission.
