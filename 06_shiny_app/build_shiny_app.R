# =============================================================================
# 06_shiny_app/build_shiny_app.R
#
# Builds the ShinyCell app
#   https://crem-bu.shinyapps.io/scRNAseq-TGFB-treatment-iAT2s/
# Fig 3D (GFP, KRT17) and Fig S3E/S3F (ASCL1, CHGA, SOX2/4/6) were exported
# from this app; 02_scRNA_analysis/fig3.R and figS3.R make equivalent plots.
#
# Sources: runShinyCell.R (CReM bioinformatics core) and the colour settings in
#          25_06_25_lauren_RNA-cluster_annotations.R (L. Ayers).
# After building, the app's server.R / ui.R were replaced with the extended
# CReM version: https://github.com/CReM-BU/cellxplorer
#
# Input : RNA_ANNOTATED_RDS
# Output: scRNAseq-TGFB-treatment-iAT2s/ (sc1conf.rds, sc1meta.rds, sc1gexpr.h5,
#         sc1gene.rds, sc1def.rds, server.R, ui.R)
#
# Run from repo root:  Rscript 06_shiny_app/build_shiny_app.R
# =============================================================================

source("config.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(ShinyCell)
  library(data.table)
})

APP_NAME <- "scRNAseq-TGFB-treatment-iAT2s"

seu <- readRDS(RNA_ANNOTATED_RDS)
DefaultAssay(seu) <- "SCT"

# Cell-type labels shown in the app (only difference from the paper: cluster 9)
seu$Cluster_annotation <- as.character(seu$subcluster_names)
seu$Cluster_annotation[seu$Cluster_annotation == "Non-lung contaminant"] <- "NKX2.1 negative"
app_cluster_cols <- CLUSTER_COLS
names(app_cluster_cols)[names(app_cluster_cols) == "Non-lung contaminant"] <- "NKX2.1 negative"

scConf <- createConfig(seu)

# Colours: matched by name to each column's level order in the config
set_cols <- function(conf, meta, cols) {
  lv <- strsplit(conf[ID == meta]$fID, "\\|")[[1]]
  modColours(conf, meta.to.mod = meta, new.colours = unname(cols[lv]))
}
scConf <- set_cols(scConf, "Cluster_annotation", app_cluster_cols)
scConf <- set_cols(scConf, "orig.ident", SAMPLE_COLS)

makeShinyApp(seu, scConf, gex.assay = "SCT", gene.mapping = TRUE,
             shiny.title = paste(
               "iPSC-derived alveolar epithelial type 2 (iAT2) cells exposed to TGFB signaling",
               "transcriptomically resemble KRT5-/KRT17+ cells observed in human lung disease."),
             shiny.footnotes = "Center for Regenerative Medicine",
             shiny.dir = APP_NAME)
