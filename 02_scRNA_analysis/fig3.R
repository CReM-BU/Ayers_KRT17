# =============================================================================
# 02_scRNA_analysis/fig3.R — Figure 3 (scRNA-seq)
#
# Panel  Output file (PLOT_DIR)                       Notes
# -----  -------------------------------------------  ---------------------------
# 3A     —                                            schematic (BioRender)
# 3B     fig3B_umap_sample.pdf
# 3C     fig3C_umap_celltype.pdf
# 3D     fig3D_umap_GFP_KRT17.pdf                     NEW script; published panel
#                                                     exported from the Shiny app (06_shiny_app/)
# 3E     fig3E_dotplot_markers.pdf
# 3F     fig3F_vln_habermann_AT2_KRT17.pdf
# 3G,3H  —                                            scTOP, see 05_scTOP/
# 3I     fig3I_gsea_hallmark_KRT17_vs_iAT2.pdf (+ .tsv in TABLE_DIR)
# 3J     fig3J_heatmap_top20_KRT17_markers.pdf
# 3K     fig3K_vln_IL32.pdf
# 3L     —                                            MISSING (Habermann et al. data;
#                                                     published panel from Shiny app)
#
# Inputs: RNA_ANNOTATED_RDS, outputs of 02_markers_and_DEG.R
# Run from repo root:  Rscript 02_scRNA_analysis/fig3.R
# =============================================================================

source("config.R")
suppressPackageStartupMessages({
  library(Seurat)      # v5.0.2
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(data.table)
  library(Matrix)
  library(pheatmap)
  library(msigdbr)
  library(fgsea)
})

rna <- readRDS(RNA_ANNOTATED_RDS)
DefaultAssay(rna) <- "SCT"

# Cluster order used for violin plots
VLN_ORDER <- c("iAT2", "Pre KRT5-/KRT17+", "KRT5-/KRT17+", "iAT1",
               "iBC 1", "iBC 2", "PNEC-like", "Proliferating (distal)",
               "SOX high TGFB-treated", "Proliferating (proximal)",
               "Differentiating Airway", "Non-lung contaminant")

# Cluster order used for heatmap columns
HEATMAP_ORDER <- c("iAT2", "Pre KRT5-/KRT17+", "KRT5-/KRT17+", "iAT1",
                   "iBC 1", "iBC 2", "Differentiating Airway",
                   "SOX high TGFB-treated", "PNEC-like",
                   "Proliferating (distal)", "Proliferating (proximal)",
                   "Non-lung contaminant")


# ---- Fig 3B: UMAP coloured by sample -----------------------------------------------
pdf(plot_path("fig3B_umap_sample.pdf"), width = 10, height = 11)
print(
  DimPlot(rna, group.by = "sample", cols = SAMPLE_COLS, pt.size = 1.5) &
    theme_umap() & theme(legend.position = "none")
)
dev.off()


# ---- Fig 3C: UMAP coloured by cell type --------------------------------------------
pdf(plot_path("fig3C_umap_celltype.pdf"), width = 10, height = 11)
print(
  DimPlot(rna, group.by = "subcluster_names", label = TRUE, cols = CLUSTER_COLS,
          repel = TRUE, label.size = 8, pt.size = 1.3) +
    theme_umap() + theme(legend.position = "none")
)
dev.off()


# ---- Fig 3D: GFP and KRT17 transcript expression -----------------------------------
# NEW: the published panel was exported from the Shiny app
# (https://crem-bu.shinyapps.io/scRNAseq-TGFB-treatment-iAT2s/). This script
# reproduces it from the same object. KRT17-GFP reads map to the EGFP_EXTRA
# transcript added in 00_alignment/scRNA/00_mkref.sh ("EGFP-EXTRA" in Seurat).
# Colour scale = ShinyCell "White-Red"; high-expressing cells plotted on top.
pdf(plot_path("fig3D_umap_GFP_KRT17.pdf"), width = 10, height = 5)
print(
  FeaturePlot(rna, features = c("EGFP-EXTRA", "KRT17"), order = TRUE, pt.size = 0.8,
              cols = SHINY_EXPR_COLS) &
    theme_umap()
)
dev.off()


# ---- Fig 3E: marker dot plot -----------------------------------------------------------
dotplot_genes <- c(
  "NKX2-1",
  "AGER", "CAV1", "EMP2",                                        # AT1
  "SFTPC", "SFTPB", "ETV5",                                      # AT2
  "ICAM1", "PTGS2", "TM4SF1", "SFN", "CLDN4", "TPM1",            # KRT5-/KRT17+
  "ITGB6", "COL1A1", "TAGLN", "FN1", "ITGAV", "SMAD3", "SOX4", "TGFB1",
  "PLAU", "CDKN2B",
  "TP53", "GADD45A", "CASP3", "TNF",
  "KRT17", "KRT7", "KRT8", "KRT19",
  "ITGB4", "S100A2", "TP63", "BCAM", "KRT5",                     # basal
  "SCGB1A1", "CCNO", "FOXJ1",                                    # differentiating airway
  "TOP2A", "MKI67", "CCNB1"                                      # proliferation
)
dotplot_rows <- c("iAT1", "iAT2", "Pre KRT5-/KRT17+", "KRT5-/KRT17+",
                  "SOX high TGFB-treated", "PNEC-like", "iBC 1", "iBC 2",
                  "Differentiating Airway", "Proliferating (proximal)",
                  "Proliferating (distal)")

Idents(rna) <- "subcluster_names"
pdf(plot_path("fig3E_dotplot_markers.pdf"), height = 6, width = 17)
print(
  DotPlot(rna, features = dotplot_genes, col.min = -1, col.max = 5, dot.scale = 10) +
    scale_y_discrete(limits = dotplot_rows) +
    RotatedAxis() + xlab("") + ylab("")
)
dev.off()


# ---- Fig 3F: Habermann AT2 and KRT5-/KRT17+ gene-set scores ------------------------
# Hb_* module scores are added in 01_scRNA_preprocessing/10x.combined.Rmd
Idents(rna) <- factor(rna$subcluster_names, levels = VLN_ORDER)
pdf(plot_path("fig3F_vln_habermann_AT2_KRT17.pdf"), width = 8, height = 5)
print(
  VlnPlot(rna, features = c("Hb_AT2", "Hb_KRT5neg.KRT17pos"),
          idents = c("iAT2", "Pre KRT5-/KRT17+", "KRT5-/KRT17+", "iAT1", "iBC 1", "iBC 2"),
          pt.size = 0, cols = CLUSTER_COLS)
)
dev.off()


# ---- Fig 3G, 3H: scTOP -- see 05_scTOP/ ---------------------------------------------


# ---- Fig 3I: Hallmark GSEA, KRT5-/KRT17+ vs iAT2 ------------------------------------
krt17.vs.iat2 <- read.table(table_path("KRT5neg_KRT17pos_vs_iAT2.markers.sct.tsv"), sep = "\t")
ranks <- setNames(krt17.vs.iat2$avg_log2FC, rownames(krt17.vs.iat2))

msig <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, gene_symbol)
pathways <- split(msig$gene_symbol, msig$gs_name)

# fgsea p-values use random permutations; seed added for reproducibility
# (the original run had no seed, so padj may differ in the last decimal places)
set.seed(1)
fgsea_res <- fgsea(pathways = pathways, stats = ranks, minSize = 5, maxSize = 500)
fgsea_res <- fgsea_res %>% arrange(dplyr::desc(abs(NES)))
fwrite(fgsea_res, file = table_path("KRT5neg_KRT17pos_vs_iAT2.hallmark.tsv"), sep = "\t")

gsea_plot <- fgsea_res %>%
  transmute(pathway = str_replace_all(pathway, "^HALLMARK_", "") |> str_replace_all("_", " "),
            NES = NES, padj = padj) %>%
  dplyr::filter(padj < 0.05) %>%            # only significant pathways are shown
  mutate(mlog10p = -log10(padj),
         pathway = reorder(pathway, NES))

pdf(plot_path("fig3I_gsea_hallmark_KRT17_vs_iAT2.pdf"), width = 8, height = 7)
print(
  ggplot(gsea_plot, aes(x = NES, y = pathway)) +
    geom_segment(aes(x = 0, xend = NES, yend = pathway), linewidth = 0.6, color = "grey80") +
    geom_point(aes(size = mlog10p, color = NES), shape = 16) +
    scale_color_gradient2(low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0,
                          limits = c(-3, 3), oob = scales::squish, name = "NES") +
    scale_size_area(max_size = 10, name = "-log10(padj)") +
    geom_vline(xintercept = 0, linetype = 3, color = "grey70") +
    labs(x = "Normalized Enrichment Score (NES)", y = NULL,
         title = "Enriched HALLMARK Pathways: \n KRT5-/KRT17+ vs. iAT2") +
    theme_minimal(base_size = 12) +
    theme(axis.text.y = element_text(size = 12), legend.position = "right")
)
dev.off()


# ---- Fig 3J: top 20 KRT5-/KRT17+ markers, mean expression per cluster -------------
# Row-scaled (z-score, capped at +/-2) mean SCT expression per cluster.
cluster_mean_heatmap <- function(obj, genes, col_order, file, width, height, cap_z = 2) {
  genes_use <- intersect(genes, rownames(obj))
  clusters  <- setNames(as.character(obj$subcluster_names), colnames(obj))
  mat <- GetAssayData(obj, assay = "SCT", layer = "data")[genes_use, , drop = FALSE]

  avg <- sapply(sort(unique(clusters)), function(cl) {
    Matrix::rowMeans(mat[, names(clusters)[clusters == cl], drop = FALSE])
  })
  avg <- matrix(avg, nrow = length(genes_use),
                dimnames = list(genes_use, sort(unique(clusters))))

  sdv <- apply(avg, 1, sd, na.rm = TRUE); sdv[sdv == 0] <- 1
  z <- (avg - rowMeans(avg, na.rm = TRUE)) / sdv
  z[z >  cap_z] <-  cap_z
  z[z < -cap_z] <- -cap_z

  pheatmap(z[, col_order], cluster_rows = TRUE, cluster_cols = FALSE,
           angle_col = 45, border_color = NA,
           filename = file, width = width, height = height)
}

top20 <- read.table(table_path("cluster_names_top20_markers.tsv"), sep = "\t")
krt17_top20 <- top20$gene[top20$cluster == "KRT5-/KRT17+"]
cluster_mean_heatmap(rna, krt17_top20, HEATMAP_ORDER,
                     file = plot_path("fig3J_heatmap_top20_KRT17_markers.pdf"),
                     width = 7, height = 5)


# ---- Fig 3K: IL32 expression per cluster -------------------------------------------
Idents(rna) <- factor(rna$subcluster_names, levels = VLN_ORDER)
pdf(plot_path("fig3K_vln_IL32.pdf"), width = 6, height = 5)
print(VlnPlot(rna, features = "IL32", pt.size = 0, cols = CLUSTER_COLS))
dev.off()


# ---- Fig 3L: IL32 in Habermann et al. epithelium -- MISSING --------------------------
# Published panel was exported from https://crem-bu.shinyapps.io/19_09_10_darrell_habermann/
# TODO: add script that plots IL32 from the Habermann et al. Seurat object.

sessionInfo()
