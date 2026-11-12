# =============================================================================
# 02_scRNA_analysis/02_markers_and_DEG.R
#
# Differential expression tables used by Fig 3I, Fig 3J and Fig S3D.
#
# Input : RNA_ANNOTATED_RDS
# Output (TABLE_DIR):
#   cluster_names_allmarkers.tsv    all cluster markers (MAST)      -> Fig S3D, Fig 3J
#   cluster_names_top20_markers.tsv top 20 per cluster (log2FC > 1) -> Fig S3D, Fig 3J
#   cluster_names_top50_markers.tsv top 50 per cluster (log2FC > 1) -> supplementary table
#   KRT5neg_KRT17pos_vs_iAT2.markers.sct.tsv  KRT5-/KRT17+ vs iAT2  -> Fig 3I (GSEA)
#
# Methods: "Cluster specific genes were calculated using MAST framework in
#           Seurat wrapper."
#
# Run from repo root:  Rscript 02_scRNA_analysis/02_markers_and_DEG.R
# =============================================================================

source("config.R")
suppressPackageStartupMessages({
  library(Seurat)   # v5.0.2
  library(MAST)
  library(dplyr)
})

rna <- readRDS(RNA_ANNOTATED_RDS)
DefaultAssay(rna) <- "SCT"
rna <- PrepSCTFindMarkers(rna)
Idents(rna) <- "subcluster_names"

# ---- Markers for every cluster ---------------------------------------------------
cluster.markers <- FindAllMarkers(rna, test.use = "MAST", only.pos = FALSE,
                                  min.pct = 0.1, verbose = FALSE)

top_n_markers <- function(markers, n) {
  markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = n) %>%
    ungroup()
}

write.table(cluster.markers,              table_path("cluster_names_allmarkers.tsv"),   sep = "\t", quote = FALSE)
write.table(top_n_markers(cluster.markers, 20), table_path("cluster_names_top20_markers.tsv"), sep = "\t", quote = FALSE)
write.table(top_n_markers(cluster.markers, 50), table_path("cluster_names_top50_markers.tsv"), sep = "\t", quote = FALSE)

# ---- KRT5-/KRT17+ vs iAT2 (input to Hallmark GSEA, Fig 3I) ----------------------
krt17.vs.iat2 <- FindMarkers(rna, ident.1 = "KRT5-/KRT17+", ident.2 = "iAT2",
                             test.use = "MAST", only.pos = FALSE,
                             min.pct = 0.1, verbose = FALSE)
write.table(krt17.vs.iat2, table_path("KRT5neg_KRT17pos_vs_iAT2.markers.sct.tsv"),
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
