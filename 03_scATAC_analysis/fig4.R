# =============================================================================
# 03_scATAC_analysis/fig4.R — Figure 4 (scATAC-seq)
#
# Panel  Output file (PLOT_DIR)              Notes
# -----  ----------------------------------  --------------------------------------
# 4A     fig4A_atac_umap_sample.pdf          NEW plotting code (not saved in original repo)
# 4B     fig4B_atac_umap_celltype.pdf        NEW plotting code (not saved in original repo)
# 4C     fig4C_atac_marker_peak_heatmap.pdf
# 4D     —                                   SCENIC+ eRegulon dot-heatmap: MISSING
#                                            plotting code; see 04_scenicplus/
# 4E-G   —                                   flow cytometry / Prism (not in repo)
#
# Input : ARCHR_PROJECT (after 02_integrate_cluster_annotate.R)
# Run from repo root:  Rscript 03_scATAC_analysis/fig4.R
# =============================================================================

source("config.R")
suppressPackageStartupMessages({
  library(ArchR)   # v1.0.3
  library(ggplot2)
})
addArchRThreads(threads = 2)

proj <- loadArchRProject(ARCHR_PROJECT)

ATAC_SAMPLE_COLS <- c(iAT1 = "#D51F26", iAT2 = "#272E6A", iBC = "#208A42", TGFB_iAT2 = "#89288F")


# ---- Fig 4A: UMAP coloured by sample (NEW) ----------------------------------------
pdf(plot_path("fig4A_atac_umap_sample.pdf"), width = 6, height = 5)
print(
  plotEmbedding(proj, colorBy = "cellColData", name = "Sample", pal = ATAC_SAMPLE_COLS) +
    theme_seurat() + xlab("umap_1") + ylab("umap_2") + ggtitle("")
)
dev.off()


# ---- Fig 4B: UMAP coloured by RNA-derived cell-type label (NEW) --------------------
pdf(plot_path("fig4B_atac_umap_celltype.pdf"), width = 6, height = 5)
print(
  plotEmbedding(proj, colorBy = "cellColData", name = "ClusterLabels",
                pal = CLUSTER_COLS[unique(proj$ClusterLabels)]) +
    theme_seurat() + xlab("umap_1") + ylab("umap_2") + ggtitle("")
)
dev.off()


# ---- Fig 4C: marker peak heatmap per cell type ----------------------------------
# Methods: FDR < 0.01 and log2FC >= 1 for marker peaks; bias-normalised for
# TSS enrichment and log10(nFrags). The heatmap display uses a looser cut-off.
markerPeaks <- getMarkerFeatures(
  ArchRProj = proj, useMatrix = "PeakMatrix", groupBy = "ClusterLabels",
  bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon"
)
markerList <- getMarkers(markerPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
saveRDS(markerPeaks, file.path(DATA_DIR, "atac_markerPeaks_ClusterLabels.rds"))

heatmapPeaks <- plotMarkerHeatmap(seMarker = markerPeaks,
                                  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
                                  transpose = TRUE)
pdf(plot_path("fig4C_atac_marker_peak_heatmap.pdf"), width = 8, height = 6)
ComplexHeatmap::draw(heatmapPeaks)
dev.off()


# ---- Fig 4D: eRegulon dot-heatmap -- MISSING (see 04_scenicplus/) ----------------

sessionInfo()
