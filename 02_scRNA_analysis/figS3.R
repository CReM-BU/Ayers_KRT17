# =============================================================================
# 02_scRNA_analysis/figS3.R — Supplementary Figure S3 (scRNA-seq)
#
# Panel  Output file (PLOT_DIR)                     Notes
# -----  -----------------------------------------  -----------------------------
# S3A    figS3A_umap_louvain_res0.25.pdf
# S3B    figS3B_umap_res0.25_subclustered.pdf
# S3C    figS3C_vln_habermann_genesets.pdf
# S3D    figS3D_heatmap_top20_DEG.pdf               page 1 = left, page 2 = right
# S3E    figS3E_umap_ASCL1_CHGA.pdf                 NEW script; published panel from Shiny app
# S3F    figS3F_umap_SOX2_SOX4_SOX6.pdf             NEW script; published panel from Shiny app
# S3G    figS3G_vln_NKX2-1.pdf
# S3H    —                                          scTOP, see 05_scTOP/
#
# Inputs: RNA_ANNOTATED_RDS, outputs of 02_markers_and_DEG.R
# Run from repo root:  Rscript 02_scRNA_analysis/figS3.R
# =============================================================================

source("config.R")
suppressPackageStartupMessages({
  library(Seurat)   # v5.0.2
  library(ggplot2)
  library(dplyr)
})

rna <- readRDS(RNA_ANNOTATED_RDS)
DefaultAssay(rna) <- "SCT"

VLN_ORDER <- c("iAT2", "Pre KRT5-/KRT17+", "KRT5-/KRT17+", "iAT1",
               "iBC 1", "iBC 2", "PNEC-like", "Proliferating (distal)",
               "SOX high TGFB-treated", "Proliferating (proximal)",
               "Differentiating Airway", "Non-lung contaminant")


# ---- Fig S3A: Louvain clusters, resolution 0.25 -----------------------------------
pdf(plot_path("figS3A_umap_louvain_res0.25.pdf"), width = 10, height = 11)
print(
  DimPlot(rna, group.by = "SCT_snn_res.0.25", label = TRUE, repel = TRUE,
          label.size = 8, pt.size = 1.3) +
    theme_umap() + theme(legend.position = "none")
)
dev.off()


# ---- Fig S3B: resolution 0.25 with sub-clustering of cluster 1 --------------------
pdf(plot_path("figS3B_umap_res0.25_subclustered.pdf"), width = 10, height = 11)
print(
  DimPlot(rna, group.by = "TGFB.subclustered", label = TRUE, repel = TRUE,
          label.size = 8, pt.size = 1.3) +
    theme_umap() + theme(legend.position = "none") +
    ggtitle("SCT_snn_res.0.25 with subclustering")
)
dev.off()


# ---- Fig S3C: Habermann et al. 20-gene set scores across all clusters --------------
# Hb_* module scores are added in 01_scRNA_preprocessing/10x.combined.Rmd
Idents(rna) <- "subcluster_names"
pdf(plot_path("figS3C_vln_habermann_genesets.pdf"), width = 12, height = 9)
print(
  VlnPlot(rna, features = c("Hb_AT2", "Hb_KRT5neg.KRT17pos", "Hb_AT1", "Hb_Basal", "Hb_MUC5Bpos"),
          pt.size = 0, cols = CLUSTER_COLS)
)
dev.off()


# ---- Fig S3D: heatmap of top 20 DEGs per cluster -----------------------------------
# Uses log-normalised RNA assay (scaled, regressing percent.mito) for display.
top20 <- read.table(table_path("cluster_names_top20_markers.tsv"), sep = "\t")

rna_hm <- rna
DefaultAssay(rna_hm) <- "RNA"
rna_hm <- NormalizeData(rna_hm)
rna_hm <- ScaleData(rna_hm, features = rownames(rna_hm),
                    vars.to.regress = "percent.mito", verbose = FALSE)
Idents(rna_hm) <- "subcluster_names"

# First 6 clusters (distal) on the left panel, last 6 (proximal/other) on the right.
# (Original code used rows 1:120 / 121:240, which assumes exactly 20 genes per cluster.)
cl_levels <- levels(rna_hm$subcluster_names)
left_genes  <- top20$gene[top20$cluster %in% cl_levels[1:6]]
right_genes <- top20$gene[top20$cluster %in% cl_levels[7:12]]

pdf(plot_path("figS3D_heatmap_top20_DEG.pdf"), width = 10, height = 15)
print(DoHeatmap(rna_hm, features = left_genes)  + scale_fill_gradientn(colors = c("blue", "white", "red")))
print(DoHeatmap(rna_hm, features = right_genes) + scale_fill_gradientn(colors = c("blue", "white", "red")))
dev.off()
rm(rna_hm)


# ---- Fig S3E: PNEC markers ASCL1 and CHGA (NEW) ------------------------------------
pdf(plot_path("figS3E_umap_ASCL1_CHGA.pdf"), width = 10, height = 5)
print(FeaturePlot(rna, features = c("ASCL1", "CHGA"), order = TRUE, cols = SHINY_EXPR_COLS) & theme_umap())
dev.off()


# ---- Fig S3F: SOX2, SOX4, SOX6 (NEW) -----------------------------------------------
pdf(plot_path("figS3F_umap_SOX2_SOX4_SOX6.pdf"), width = 15, height = 5)
print(FeaturePlot(rna, features = c("SOX2", "SOX4", "SOX6"), order = TRUE, ncol = 3,
                  cols = SHINY_EXPR_COLS) & theme_umap())
dev.off()


# ---- Fig S3G: NKX2-1 expression per cluster ----------------------------------------
Idents(rna) <- factor(rna$subcluster_names, levels = VLN_ORDER)
pdf(plot_path("figS3G_vln_NKX2-1.pdf"), width = 12, height = 9)
print(VlnPlot(rna, features = "NKX2-1", pt.size = 0, cols = CLUSTER_COLS))
dev.off()


# ---- Fig S3H: scTOP -- see 05_scTOP/ -------------------------------------------------

sessionInfo()
