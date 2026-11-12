# =============================================================================
# 03_scATAC_analysis/figS4.R — Supplementary Figure S4 (scATAC-seq)
#
# Panel letters follow the figure PDF (krt17_suppfigs_v8). The v9.1 legend
# letters these A, B, B, C, D, E -- the legend needs updating to A-F.
#
# Panel  Output file (PLOT_DIR)                         Notes
# -----  ---------------------------------------------  -------------------------
# S4A    figS4A_atac_umap_clusters.pdf
# S4B    figS4B_label_transfer_heatmap_clusters.pdf     before sub-clustering
# S4C    figS4C_atac_umap_subclusters.pdf               after sub-clustering of C7
# S4D    figS4D_label_transfer_heatmap_subclusters.pdf  with C7 sub-clusters
# S4E    figS4E_<TF>_deviation_ridge.pdf                ridge plots, 6 motifs
#        figS4E_<TF>_motif_logo.pdf                     sequence logos
#        figS4E_<TF>_deviation_umap.pdf                 UMAPs
# S4F    —                                              full eRegulon heatmap:
#                                                       MISSING plotting code
#
# Input : ARCHR_PROJECT (after 02_integrate_cluster_annotate.R)
# Run from repo root:  Rscript 03_scATAC_analysis/figS4.R
# =============================================================================

source("config.R")
suppressPackageStartupMessages({
  library(ArchR)    # v1.0.3
  library(ggplot2)
  library(pheatmap)
  library(ggseqlogo)
  library(TFBSTools)
})
addArchRThreads(threads = 2)

proj <- loadArchRProject(ARCHR_PROJECT)

# Row-normalised confusion matrix: ATAC cluster (rows) x transferred RNA label (cols)
label_heatmap <- function(clusters, file) {
  cM <- confusionMatrix(paste0(clusters), paste0(proj$predictedGroup))
  cM <- cM / Matrix::rowSums(cM)
  p <- pheatmap::pheatmap(mat = as.matrix(cM), color = paletteContinuous("whiteBlue"),
                          border_color = "black", silent = TRUE)
  pdf(file, width = 6, height = 4); grid::grid.draw(p$gtable); dev.off()
}

atac_umap <- function(name, file) {
  pdf(file, width = 6, height = 5)
  print(plotEmbedding(proj, name = name, discreteSet = "calm") +
          theme_seurat() + xlab("umap_1") + ylab("umap_2") + ggtitle(""))
  dev.off()
}


# ---- Fig S4A: Louvain clusters (resolution 0.25) ---------------------------------
atac_umap("Clusters_0.25", plot_path("figS4A_atac_umap_clusters.pdf"))

# ---- Fig S4B: label-transfer heatmap, before sub-clustering ------------------------
label_heatmap(proj$Clusters_0.25, plot_path("figS4B_label_transfer_heatmap_clusters.pdf"))

# ---- Fig S4C: clusters after sub-clustering of C7 ----------------------------------
atac_umap("subclusters", plot_path("figS4C_atac_umap_subclusters.pdf"))

# ---- Fig S4D: label-transfer heatmap, with C7 sub-clusters -------------------------
label_heatmap(proj$subclusters, plot_path("figS4D_label_transfer_heatmap_subclusters.pdf"))


# ---- Fig S4E: chromVAR motif deviations + sequence logos ----------------------------
# Motif names in the cisBP MotifMatrix
S4D_MOTIFS <- c(NKX2_1 = "z:NKX21_476", TEAD4 = "z:TEAD4_797", SMAD3 = "z:SMAD3_743",
                TP63   = "z:TP63_704",  NFKB1 = "z:NFKB1_719", FOSL1 = "z:FOSL1_142")
stopifnot(all(S4D_MOTIFS %in% getFeatures(proj, useMatrix = "MotifMatrix")))

ridge_groups <- c("Pre KRT5-/KRT17+", "KRT5-/KRT17+", "iAT1", "iAT2", "iBC 1", "iBC 2",
                  "Differentiating Airway")
proj_sub <- proj[proj$ClusterLabels %in% ridge_groups, ]
ridges <- plotGroups(ArchRProj = proj_sub, groupBy = "ClusterLabels",
                     colorBy = "MotifMatrix", name = unname(S4D_MOTIFS),
                     pal = CLUSTER_COLS)

# cisBP PWM (log-odds) -> probability matrix for ggseqlogo
PWMatrixToProbMatrix <- function(x) {
  if (class(x) != "PWMatrix") stop("x must be a TFBSTools::PWMatrix object")
  m <- exp(as(x, "matrix")) * TFBSTools::bg(x) / sum(TFBSTools::bg(x))
  t(t(m) / colSums(m))
}
motif_pwms <- getPeakAnnotation(proj, "Motif")$motifs

for (tf in names(S4D_MOTIFS)) {
  motif <- S4D_MOTIFS[[tf]]

  pdf(plot_path(paste0("figS4E_", tf, "_deviation_ridge.pdf")), height = 3, width = 3)
  print(ridges[[motif]])
  dev.off()

  pdf(plot_path(paste0("figS4E_", tf, "_motif_logo.pdf")), height = 1.5, width = 3)
  print(ggseqlogo(PWMatrixToProbMatrix(motif_pwms[[sub("^z:", "", motif)]]), method = "bits"))
  dev.off()

  pdf(plot_path(paste0("figS4E_", tf, "_deviation_umap.pdf")), height = 6, width = 5)
  print(plotEmbedding(proj, colorBy = "MotifMatrix", name = motif) +
          theme_seurat() + xlab("umap_1") + ylab("umap_2") + ggtitle(""))
  dev.off()
}


# ---- Fig S4F: full eRegulon heatmap -- MISSING (see 04_scenicplus/) ------------------

sessionInfo()
