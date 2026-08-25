# CYL/ZCP Scanpy cluster validation and annotation correction

Read this reference together with `scanpy-transcriptome.md` whenever reviewing markers, changing a cluster label, updating the annotation dotplot, or transferring notebook conclusions into the production Python script.

## Non-negotiable distinction

- `leiden` is the cluster-level validation key. Preserve one evidence row per cluster even when several clusters share a biological label.
- `cell_type` is the reviewed, merged presentation level. Marker dotplots and final cell-type UMAPs group by this field, not by `leiden`.
- In the current review, clusters 5/9 merge as `Macrophages`, 7/8 merge as `Endothelial cells`, and 16/17 merge as `NA`. Their cluster-specific evidence must not be discarded.
- Cluster numbers are run-specific. Never copy the reviewed 0–17 mapping unless the fresh observed cluster set matches exactly.

## Evidence hierarchy for annotation

Use evidence in this order:

1. A coherent multi-gene positive program in normalized expression, supported by both detection fraction and mean expression.
2. Expected exclusion or negative evidence against competing lineages; inspect raw/log-normalized expression rather than scaled values alone.
3. Cluster QC, Scrublet scores, cell number, and sample composition. Sample skew is a warning, not proof of artifact.
4. UMAP neighbourhood and separation as supporting context only.
5. Prior broad CYL/ZCP labels as lineage-level auxiliary evidence only; barcode-resolved prior labels are unavailable.

Do not auto-correct a label from one top marker, UMAP position, cluster number, or the prior broad annotation alone. A state label such as `Cycling cells` or `MT-high AT2-like` also does not by itself define a distinct lineage or justify cell removal.

## Reviewed 0–17 decisions

| Leiden | Current label | Confidence | Validation and future correction rule |
|---|---|---|---|
| 0 | AT2 | high | SFTPC/SFTPB/ABCA3/LPCAT1 form a coherent AT2 program. |
| 1 | Secretory epithelial | medium | Secretory/epithelial localization with NEDD4L/SFTA3 and retained surfactant signal; do not force a narrower subtype without a coherent secretory panel. |
| 2 | Fibroblasts | high | COL1A2/COL3A1/COL5A1/COL6A3/PDGFRA support fibroblast identity. |
| 3 | Ciliated cells | high | DNAH/CFAP/HYDIN/FOXJ1-family ciliogenesis program. |
| 4 | Secretory / mucous epithelial | high | BPIFB1/MUC4/ERN2/TMC5 support secretory-mucous identity. |
| 5 | Macrophages | high | CD163/MRC1/CTSB/FCER1G; retain as a separate macrophage state at cluster level. |
| 6 | AT1-like | medium | CAV1/HOPX/CAV2 are about 65/63/36%, but AGER/PDPN/AQP5 only about 15/11/3%, while SFTPB/LPCAT1/ABCA3/SFTPC remain about 96/73/49/33%. Do not promote to pure AT1 until the canonical AT1 program strengthens and AT2 carry-over decreases. |
| 7 | Endothelial cells | high | EPAS1/PECAM1/VWF/BTNL9; retain its endothelial state separately from cluster 8 in validation tables. |
| 8 | Endothelial cells | high | VWF/PTPRB/PECAM1/EPAS1; merge with 7 only at primary cell-type level. |
| 9 | Macrophages | high | PPARG/MRC1/CD163/MSR1 suggest a macrophage state; merge with 5 for primary annotation, then subset if subtype analysis is requested. |
| 10 | Basal cells | high | EGFR/KRT15/COL7A1/TP63/KRT5 provide a coherent basal program. |
| 11 | Smooth muscle / mural cells | high | MYH11/LMOD1/CARMN/PDGFRB/COL4A1 support mural identity. |
| 12 | MT-high AT2-like | medium | Mitochondrial-high state retains SFTPC/SFTPA2/SFTPB/ABCA3. Keep as a QC/state label; do not automatically delete it. |
| 13 | T cells | high | PTPRC/CD247/BCL11B/ITK/DOCK2 support T lineage. Prior broad `NA` cannot override this program. |
| 14 | Cycling cells | medium | DIAPH3/FANCI/MELK/RRM2/ECT2/ANLN support cycling. It lies near epithelial cells, but only call `Cycling epithelial cells` after demonstrating EPCAM/KRT8/KRT18/KRT19 coherently. |
| 15 | Lymphatic endothelial cells | high | PROX1/FLT4/LYVE1/CCL21/MMRN1/RELN are detected in about 29/40/22/19/59/49%, with DCN 0%; pan-endothelial context supports lymphatic endothelium. Prior `NA` cannot override this evidence. |
| 16 | NA | low | 50 cells, about 94% CYL; epithelial-like clues but no complete lineage or cycling program. Keep `NA` until independent markers cohere. |
| 17 | NA | low | 20 cells, about 90% ZCP; stromal-like COL1A2/DCN/COL3A1 clues but too small/sample-skewed for a final label. Keep `NA`. |

The prior broad annotation is documented in `Supplementary/prior_annotation_crosswalk.tsv`. It supports the major lineage structure, but it is not a cell-by-cell confusion matrix.

## Required correction workflow

1. Before execution, validate the internally reviewed mapping and marker panel:

   ```bash
   /home/lijia/jiangyuanpei/miniforge3/envs/allcools/bin/python \
     Scripts/Scanpy/Scripts/02_validate_annotation_config.py
   ```

2. Run only into a new empty `Results/Scanpy/<run_name>/` root. Inspect `annotation_guard_status.json` before treating `cell_type` as reviewed.
3. If observed IDs differ from 0–17, do not renumber or nearest-neighbour-transfer the old mapping. All cells must remain `Unassigned`; use `tables/leiden_ranked_markers.tsv.gz`, `tables/cluster_annotation_validation.tsv`, the ranked-marker figure, cluster QC/sample fractions, and lineage-marker evidence to review every new cluster.
4. For an uncertain cluster, compare at least the proposed lineage, its closest competing lineage, QC/doublet metrics, and sample composition. Use `NA` or `Unassigned` when evidence is incomplete.
5. Update the notebook interpretation first, then synchronize `CLUSTER_ANNOTATIONS`, `DOTPLOT_MARKERS`, `CELL_TYPE_ORDER`, `Scripts/Scanpy/Report.md`, and this reference. Run the validator again.
6. Keep three to four non-duplicated markers per merged cell type. The annotation dotplot must group by `cell_type`, and its y-axis labels and top marker-group labels must use the same dendrogram-synchronized order.

## Machine-readable evidence expected from production

- `annotation_guard_status.json`: exact observed-versus-reviewed cluster-set comparison.
- `tables/reviewed_cluster_annotations.tsv`: reviewed label, confidence, and concise evidence per cluster.
- `tables/cluster_annotation_validation.tsv`: cluster annotation joined to QC and sample composition for the current run.
- `tables/cell_type_cluster_membership.tsv`: explicit mapping from merged cell type back to member Leiden clusters.
- `tables/dotplot_marker_panel.tsv`: ordered marker groups and gene-presence status.
- `tables/leiden_ranked_markers.tsv.gz`: full ranked-marker evidence used for re-review.

The configuration validator checks structural consistency only. A passing result does not validate biological identity or prove that a fresh run reproduced the notebook result.
