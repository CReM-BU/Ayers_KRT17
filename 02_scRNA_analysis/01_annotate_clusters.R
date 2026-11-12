# =============================================================================
# 02_scRNA_analysis/01_annotate_clusters.R
#
# Sub-clusters cluster 1 (TGFB-treated cells) and assigns the cell-type names
# used throughout the paper.
#
# Source: 25_06_25_lauren_RNA-cluster_annotations.R (L. Ayers)
#
# Input : RNA_RAW_RDS        (combined Seurat object from 01_scRNA_preprocessing:
#                             combined/rds/sc.Rds, shared as RNA.rds)
# Output: RNA_ANNOTATED_RDS  (adds `TGFB.subclustered` and `subcluster_names`)
#
# Clusters are from SCT_snn_res.0.25 (Louvain, resolution 0.25; Fig S3A).
# Methods: "Subclustering was performed on cluster 1 using Seurat's
#           FindSubCluster() function (resolution = 0.25)."
#
# Run from repo root:  Rscript 02_scRNA_analysis/01_annotate_clusters.R
# =============================================================================

source("config.R")
suppressPackageStartupMessages({
  library(Seurat)   # v5.0.2
  library(ggplot2)
})

rna <- readRDS(RNA_RAW_RDS)
DefaultAssay(rna) <- "SCT"

# ---- Sub-cluster cluster 1 at resolution 0.25 ---------------------------------
Idents(rna) <- "SCT_snn_res.0.25"
rna <- FindSubCluster(rna, cluster = "1", subcluster.name = "TGFB.subclustered",
                      graph.name = "SCT_nn", resolution = 0.25, algorithm = 1)

# Sub-cluster 1_3 contains only 2 cells; merge it into 1_1
rna$TGFB.subclustered <- ifelse(rna$TGFB.subclustered == "1_3", "1_1", rna$TGFB.subclustered)

# ---- Name clusters --------------------------------------------------------------
Idents(rna) <- "TGFB.subclustered"
rna <- RenameIdents(
  rna,
  "0"   = "iAT1",
  "1_0" = "SOX high TGFB-treated",
  "1_1" = "KRT5-/KRT17+",
  "1_2" = "Pre KRT5-/KRT17+",
  "2"   = "Proliferating (distal)",
  "3"   = "iAT2",
  "4"   = "iBC 1",
  "5"   = "iBC 2",
  "6"   = "Differentiating Airway",
  "7"   = "Proliferating (proximal)",
  "8"   = "PNEC-like",
  "9"   = "Non-lung contaminant"   # "NKX2-1 negative" in the v8 figure PDFs and the Shiny app
)
rna$subcluster_names <- Idents(rna)

stopifnot(all(levels(rna$subcluster_names) %in% names(CLUSTER_COLS)))
print(table(rna$subcluster_names, rna$sample))

saveRDS(rna, RNA_ANNOTATED_RDS)
