# =============================================================================
# 04_scenicplus/01_export_inputs.R
#
# Exports scATAC (ArchR) and scRNA (Seurat) data in the formats needed by
# pycisTopic / SCENIC+.
#
# Inputs : ARCHR_PROJECT, RNA_ANNOTATED_RDS
# Outputs (SCENIC_DIR):
#   atac/peaks.consensus.bed                       ArchR consensus peak set
#   atac/pycisTopic_input/peak_counts.mtx          peak x cell counts
#   atac/pycisTopic_input/{peaks,cells}.tsv
#   atac/pycisTopic_input/cell_metadata.csv
#   rna/rna.h5ad                                   RNA counts + log-normalised data
#   markerpeaks/<subcluster>.bed, DApeaks_union.bed  DARs per sub-cluster (not used by
#                                                  the final SCENIC+ run, which uses topic region sets)
#   markerpeaks/ArchR_allPeaks.bed                 all peaks -> cistarget database
#
# Run from repo root:  Rscript 04_scenicplus/01_export_inputs.R
# =============================================================================

source("config.R")
suppressPackageStartupMessages({
  library(ArchR)                 # v1.0.3
  library(Seurat)                # v5.0.2
  library(SingleCellExperiment)
  library(zellkonverter)
  library(data.table)
  library(dplyr)
})

atac_in  <- file.path(SCENIC_DIR, "atac", "pycisTopic_input")
rna_out  <- file.path(SCENIC_DIR, "rna")
peak_out <- file.path(SCENIC_DIR, "markerpeaks")
for (d in c(atac_in, rna_out, peak_out)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

proj <- loadArchRProject(ARCHR_PROJECT)

# ---- ATAC: consensus peaks and peak x cell matrix ---------------------------------
peaks_gr <- getPeakSet(proj)
peaks_df <- data.frame(seqnames = as.character(seqnames(peaks_gr)),
                       start = start(peaks_gr), end = end(peaks_gr))
peaks_df$name <- paste0(peaks_df$seqnames, ":", peaks_df$start, "-", peaks_df$end)
write.table(peaks_df[, c("seqnames", "start", "end", "name")],
            file.path(SCENIC_DIR, "atac", "peaks.consensus.bed"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

pm  <- getMatrixFromProject(proj, useMatrix = "PeakMatrix")
mat <- assays(pm)[["PeakMatrix"]]   # dgCMatrix, peaks x cells
rownames(mat) <- peaks_df$name      # "chr:start-end", same order as getPeakSet()
Matrix::writeMM(mat, file.path(atac_in, "peak_counts.mtx"))
write.table(rownames(mat), file.path(atac_in, "peaks.tsv"), quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(colnames(mat), file.path(atac_in, "cells.tsv"), quote = FALSE, row.names = FALSE, col.names = FALSE)

cellmeta <- as.data.frame(getCellColData(proj))
cellmeta$cell <- rownames(cellmeta)
write.csv(cellmeta, file.path(atac_in, "cell_metadata.csv"), row.names = FALSE)

# ---- ATAC: differentially accessible regions per sub-cluster ----------------------
markersPeaks <- getMarkerFeatures(ArchRProj = proj, useMatrix = "PeakMatrix",
                                  groupBy = "subclusters", testMethod = "wilcoxon",
                                  bias = c("TSSEnrichment", "log10(nFrags)"))
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 1")

all_beds <- list()
for (cl in names(markerList)) {
  dt  <- as.data.table(markerList[[cl]])
  bed <- unique(data.table(seqnames = dt$seqnames, start = dt$start, end = dt$end)[
    !is.na(seqnames) & !is.na(start) & !is.na(end)])
  all_beds[[cl]] <- bed
  fwrite(bed, file.path(peak_out, paste0(cl, ".bed")), sep = "\t", col.names = FALSE)
}
union_bed <- unique(rbindlist(all_beds, fill = TRUE))
setorder(union_bed, seqnames, start, end)
fwrite(union_bed, file.path(peak_out, "DApeaks_union.bed"), sep = "\t", col.names = FALSE)

# All peaks -> input to 02_prepare_padded_fasta.sh (cistarget database)
fwrite(peaks_df[, c("seqnames", "start", "end")], file.path(peak_out, "ArchR_allPeaks.bed"),
       sep = "\t", col.names = FALSE)

# ---- RNA: Seurat -> AnnData (.h5ad) ------------------------------------------------
# NOTE: original script read 'rna-tgfb-subclusters-100925.rds'; now uses the
# annotated object. Idents are set to cell type because SCENIC+ groups cells by
# GEX:ident (05_scenicplus_prepare_data.sh).
rna <- readRDS(RNA_ANNOTATED_RDS)
DefaultAssay(rna) <- "RNA"
rna <- JoinLayers(rna)
rna <- NormalizeData(rna, assay = "RNA")
Idents(rna) <- "subcluster_names"

sce <- as.SingleCellExperiment(rna, assay = "RNA")
assay(sce, "counts")    <- GetAssayData(rna, assay = "RNA", layer = "counts")
assay(sce, "logcounts") <- GetAssayData(rna, assay = "RNA", layer = "data")
reducedDim(sce, "X_pca")  <- Embeddings(rna, "pca")
reducedDim(sce, "X_umap") <- Embeddings(rna, "umap")

zellkonverter::writeH5AD(sce, file = file.path(rna_out, "rna.h5ad"), X_name = "counts")
