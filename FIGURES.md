# Figure panel → code map

Only computational panels are covered. Bar plots (GraphPad Prism), flow cytometry
plots, microscopy and schematics are not part of this repository.

**Status:** ✅ code in repo · 🆕 new plotting code added during clean-up (please run and
check it matches the published panel) · ❌ missing · ⏳ pending from collaborator

Output files are written to `PLOT_DIR` (see `config.R`).

Panel letters follow the figure PDFs (`krt17_mainfigs_v8`, `krt17_suppfigs_v8`).
Panels marked "Shiny" were exported from the ShinyCell app (`06_shiny_app/`); the
scripts listed here make an equivalent plot from the same object.

## Figure 3 — scRNA-seq

| Panel | Description | Script | Output | Status |
|---|---|---|---|---|
| 3B | UMAP by sample | `02_scRNA_analysis/fig3.R` | `fig3B_umap_sample.pdf` | ✅ |
| 3C | UMAP by cell type | `02_scRNA_analysis/fig3.R` | `fig3C_umap_celltype.pdf` | ✅ |
| 3D | GFP and KRT17 expression | `02_scRNA_analysis/fig3.R` | `fig3D_umap_GFP_KRT17.pdf` | 🆕 (Shiny) |
| 3E | Marker dot plot | `02_scRNA_analysis/fig3.R` | `fig3E_dotplot_markers.pdf` | ✅ |
| 3F | Habermann AT2 / KRT5-/KRT17+ gene-set scores | `02_scRNA_analysis/fig3.R` | `fig3F_vln_habermann_AT2_KRT17.pdf` | ✅ |
| 3G, 3H | scTOP alignment scores | `05_scTOP/` | — | ⏳ |
| 3I | Hallmark GSEA, KRT5-/KRT17+ vs iAT2 | `02_markers_and_DEG.R` → `fig3.R` | `fig3I_gsea_hallmark_KRT17_vs_iAT2.pdf` | ✅ |
| 3J | Top 20 KRT5-/KRT17+ markers across clusters | `02_markers_and_DEG.R` → `fig3.R` | `fig3J_heatmap_top20_KRT17_markers.pdf` | ✅ |
| 3K | IL32 across clusters | `02_scRNA_analysis/fig3.R` | `fig3K_vln_IL32.pdf` | ✅ |
| 3L | IL32 in Habermann et al. epithelium | — | — | ❌ (Habermann Shiny app) |

## Figure 4 — scATAC-seq

| Panel | Description | Script | Output | Status |
|---|---|---|---|---|
| 4A | UMAP by sample | `03_scATAC_analysis/fig4.R` | `fig4A_atac_umap_sample.pdf` | 🆕 |
| 4B | UMAP by RNA-derived label (with C7 sub-clusters) | `03_scATAC_analysis/fig4.R` | `fig4B_atac_umap_celltype.pdf` | 🆕 |
| 4C | Marker peak heatmap | `03_scATAC_analysis/fig4.R` | `fig4C_atac_marker_peak_heatmap.pdf` | ✅ |
| 4D | Cell-type eRegulons (SCENIC+) | `04_scenicplus/` | `04_scenicplus/results/eregulon_*_heatmap.pdf` (outputs only) | ❌ plotting code + later SCENIC+ steps |

## Figure S3 — scRNA-seq

| Panel | Description | Script | Output | Status |
|---|---|---|---|---|
| S3A | Louvain clusters, res 0.25 | `02_scRNA_analysis/figS3.R` | `figS3A_umap_louvain_res0.25.pdf` | ✅ |
| S3B | Res 0.25 with cluster 1 sub-clustered | `02_scRNA_analysis/figS3.R` | `figS3B_umap_res0.25_subclustered.pdf` | ✅ |
| S3C | Five Habermann gene-set scores | `02_scRNA_analysis/figS3.R` | `figS3C_vln_habermann_genesets.pdf` | ✅ |
| S3D | Top 20 DEGs per cluster heatmap | `02_markers_and_DEG.R` → `figS3.R` | `figS3D_heatmap_top20_DEG.pdf` | ✅ |
| S3E | ASCL1, CHGA | `02_scRNA_analysis/figS3.R` | `figS3E_umap_ASCL1_CHGA.pdf` | 🆕 (Shiny) |
| S3F | SOX2, SOX4, SOX6 | `02_scRNA_analysis/figS3.R` | `figS3F_umap_SOX2_SOX4_SOX6.pdf` | 🆕 (Shiny) |
| S3G | NKX2-1 across clusters | `02_scRNA_analysis/figS3.R` | `figS3G_vln_NKX2-1.pdf` | ✅ |
| S3H | scTOP AT2 vs basal | `05_scTOP/` | — | ⏳ |

## Figure S4 — scATAC-seq

The v9.1 legend letters these panels A, B, B, C, D, E; the PDF letters them A–F (used here).

| Panel | Description | Script | Output | Status |
|---|---|---|---|---|
| S4A | Louvain clusters, res 0.25 | `03_scATAC_analysis/figS4.R` | `figS4A_atac_umap_clusters.pdf` | ✅ |
| S4B | Label-transfer heatmap, before sub-clustering | `03_scATAC_analysis/figS4.R` | `figS4B_label_transfer_heatmap_clusters.pdf` | ✅ |
| S4C | Clusters after C7 sub-clustering | `03_scATAC_analysis/figS4.R` | `figS4C_atac_umap_subclusters.pdf` | ✅ |
| S4D | Label-transfer heatmap, with C7 sub-clusters | `03_scATAC_analysis/figS4.R` | `figS4D_label_transfer_heatmap_subclusters.pdf` | ✅ |
| S4E | Motif deviations (ridge, logo, UMAP), 6 TFs | `03_scATAC_analysis/figS4.R` | `figS4E_<TF>_{deviation_ridge,motif_logo,deviation_umap}.pdf` | ✅ |
| S4F | Full eRegulon heatmap | `04_scenicplus/` | `04_scenicplus/results/eregulon_*_heatmap.pdf` (outputs only) | ❌ plotting code |

## Supplementary Figures S1, S2, S5, S6

No computational panels (gels, karyotype, flow cytometry, RT-qPCR, imaging).
