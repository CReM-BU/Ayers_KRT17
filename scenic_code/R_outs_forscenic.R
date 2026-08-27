setwd('/restricted/projectnb/crem-trainees/Kotton_Lab/Lauren/25_06_24_lauren_RDS/')

###### ATAC
# export ATAC data
library(ArchR)
proj <- loadArchRProject("LJA_250624_ATAC")

peaks_gr <- getPeakSet(proj)                   
peaks_df <- data.frame(seqnames=as.character(seqnames(peaks_gr)),
                       start=start(peaks_gr),
                       end=end(peaks_gr))
peaks_df$name <- paste0(peaks_df$seqnames, ":", peaks_df$start, "-", peaks_df$end)
write.table(peaks_df[,c("seqnames","start","end","name")],
            "../scenic/atac/peaks.consensus.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

pm <- getMatrixFromProject(proj, useMatrix="PeakMatrix")
mat <- assays(pm)[["PeakMatrix"]]   # dgCMatrix (features × cells)
rowData(pm)$name <- paste0(peaks_df$seqnames,":",peaks_df$start,"-",peaks_df$end)
peaks <- rowData(pm)$name           # should match "chr:start-end"
cells <- colnames(mat)
rownames(mat) <- peaks
dir.create("../scenic/atac/pycisTopic_input", showWarnings = FALSE)
Matrix::writeMM(mat, "../scenic/atac/pycisTopic_input/peak_counts.mtx")
write.table
write.table(peaks, "../scenic/atac/pycisTopic_input/peaks.tsv", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(cells, "../scenic/atac/pycisTopic_input/cells.tsv", quote=FALSE, row.names=FALSE, col.names=FALSE)

cellmeta <- as.data.frame(getCellColData(proj))
cellmeta$cell <- rownames(cellmeta)
write.csv(cellmeta, "../scenic/atac/pycisTopic_input/cell_metadata.csv", row.names=FALSE)

###### RNA
# make RNA seurat object into anndata
library(Seurat, lib.loc = "/share/pkg.8/r/4.4.0/install/lib64/R/library")
library(SeuratDisk)

rna <- readRDS('rna-tgfb-subclusters-100925.rds')
DefaultAssay(rna) <- "RNA"
library(SingleCellExperiment)
library(zellkonverter)

DefaultAssay(rna) <- "RNA"
rna <- NormalizeData(rna, assay = "RNA")

sce <- as.SingleCellExperiment(rna, assay = "RNA")

# map Seurat layers into assays
assay(sce, "counts")    <- GetAssayData(rna, assay = "RNA", layer = "counts")
assay(sce, "logcounts") <- GetAssayData(rna, assay = "RNA", layer = "data")

reducedDim(sce, "X_pca")  <- Embeddings(rna, "pca")
reducedDim(sce, "X_umap") <- Embeddings(rna, "umap")

# export to anndata (.h5ad)
zellkonverter::writeH5AD(
  sce,
  file = "/restricted/projectnb/crem-trainees/Kotton_Lab/Lauren/scenic/rna/rna.h5ad",
  X_name = "counts"
  )
)


### export DARs

markersPeaks <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  groupBy   = "subclusters",            
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment","log10(nFrags)")
)

markerList <- getMarkers(
  markersPeaks,
  cutOff = "FDR <= 0.05 & Log2FC >= 1"
)

outdir <- '/restricted/projectnb/crem-trainees/Kotton_Lab/Lauren/25_06_24_lauren_RDS/markerpeaks'
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

all_beds <- list()

# write each DAR peakset to bed
for (cl in names(markerList)) {
  dt <- as.data.table(markerList[[cl]])
  bed <- data.table(seqnames = dt$seqnames, start = dt$start, end = dt$end)
  bed <- unique(bed[!is.na(seqnames) & !is.na(start) & !is.na(end)])
  all_beds[[cl]] <- bed
  fwrite(bed, file.path(outdir, paste0(cl, ".bed")), sep = "\t", col.names = FALSE)
}

# make union DAR peakset
union_bed <- unique(rbindlist(all_beds, fill = TRUE))
setorder(union_bed, seqnames, start, end)
fwrite(union_bed, file.path(outdir, "DApeaks_union.bed"), sep = "\t", col.names = FALSE)


peakNames <- getFeatures(proj, useMatrix = "PeakMatrix")

peak_dt <- data.frame(seqnames = peakNames@seqnames, ranges = peakNames@ranges)

peak_dt <- peak_dt %>% dplyr::select(seqnames, ranges.start, ranges.end)

out_path <- "/restricted/projectnb/crem-trainees/Kotton_Lab/Lauren/25_06_24_lauren_RDS/markerpeaks/ArchR_allPeaks.bed"
fwrite(peak_dt, out_path, sep = "\t", col.names = FALSE)



