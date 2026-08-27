# IPF tissue MethylVI report

## Cell inventory

Input data: `Data/30wcov` (Bismark coverage files). The current MethylVI run retains all 6,554 cells.

| Cohort | Raw `.cov` cells | Cells in MethylVI |
|---|---:|---:|
| CYL | 3,165 | 3,165 |
| ZCP | 3,389 | 3,389 |
| Total | 6,554 | 6,554 |

## Manual cell-type annotation

Annotations are read from `Supplementary/manual_celltype_annotation.tsv` and are summarized below.

| Cell type | CYL | ZCP | Total |
|---|---:|---:|---:|
| AT1 | 692 | 520 | 1,212 |
| AT2 | 635 | 289 | 924 |
| Basal | 122 | 100 | 222 |
| Ciliated | 189 | 372 | 561 |
| Endothelial | 312 | 141 | 453 |
| Fibroblast | 237 | 235 | 472 |
| Goblet/Secretory | 160 | 254 | 414 |
| Macrophage/Myeloid | 309 | 363 | 672 |
| Proliferating | 73 | 30 | 103 |
| Smooth muscle/Pericyte | 58 | 20 | 78 |
| Unknown (unannotated) | 378 | 1,065 | 1,443 |
| Total | 3,165 | 3,389 | 6,554 |

## Interpretation note

ZCP contains substantially more unannotated cells (1,065/3,389) than CYL (378/3,165). This imbalance must be considered when comparing cohort composition or interpreting supervised UMAPs based on `manual_celltype`.
