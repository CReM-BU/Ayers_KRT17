# =============================================================================
# 03_scATAC_analysis/01_build_archr_project.R
#
# Builds the ArchR project from cellranger-atac fragment files: QC filtering,
# doublet removal, iterative LSI, UMAP, peak calling, motif annotation and
# chromVAR motif deviations.
#
# Input : fragments.tsv.gz from 00_alignment/scATAC/02_cellranger_atac_count.sh
# Output: ARCHR_PROJECT (ArchR project directory)
#
# Methods: TSS enrichment >= 4 and >= 1000 fragments (createArrowFiles defaults),
#          filterDoublets, iterative LSI, Louvain clustering, UMAP, MACS2 peaks
#          per cluster + iterative overlap merging, cisBP motifs, chromVAR.
#
# This code was commented out at the top of the original ATAC_fig4.R; it is
# unchanged except that addDeviationsMatrix() was added (the MotifMatrix used
# for Fig S4D is created by it, but the call was missing).
#
# Run from repo root:  Rscript 03_scATAC_analysis/01_build_archr_project.R
# =============================================================================

source("config.R")
suppressPackageStartupMessages({
  library(ArchR)   # v1.0.3
})
set.seed(1)

FRAG_DIR <- "/restricted/projectnb/crem-bioinfo/project_workspace/25_06_24_lauren/ATAC/calculations/cellranger_count"
input_files <- c(
  iAT2      = file.path(FRAG_DIR, "AL_KD_ATAC_1/outs/fragments.tsv.gz"),
  iAT1      = file.path(FRAG_DIR, "AL_KD_ATAC_2/outs/fragments.tsv.gz"),
  TGFB_iAT2 = file.path(FRAG_DIR, "AL_KD_ATAC_3/outs/fragments.tsv.gz"),
  iBC       = file.path(FRAG_DIR, "AL_KD_ATAC_4/outs/fragments.tsv.gz")
)

addArchRGenome("hg38")
addArchRThreads(threads = 1)

# Arrow files are written to the working directory; keep them with the data
setwd(DATA_DIR)

# ---- QC + doublets -------------------------------------------------------------
ArrowFiles <- createArrowFiles(inputFiles = input_files,
                               sampleNames = names(input_files))   # minTSS = 4, minFrags = 1000 (defaults)
doubScores <- addDoubletScores(input = ArrowFiles, k = 10, knnMethod = "UMAP", LSIMethod = 1)

proj <- ArchRProject(ArrowFiles = ArrowFiles, outputDirectory = basename(ARCHR_PROJECT),
                     copyArrows = TRUE)
proj <- filterDoublets(proj)

# ---- Dimensionality reduction, initial clusters, UMAP ---------------------------
proj <- addIterativeLSI(proj, useMatrix = "TileMatrix", name = "IterativeLSI")
proj <- addClusters(proj, reducedDims = "IterativeLSI")
proj <- addUMAP(proj, reducedDims = "IterativeLSI")
proj <- saveArchRProject(proj)

# ---- Peak calling ---------------------------------------------------------------
proj <- addGroupCoverages(proj)
proj <- addReproduciblePeakSet(proj, groupBy = "Clusters", pathToMacs2 = findMacs2())
proj <- addPeakMatrix(proj)
proj <- saveArchRProject(proj)

# ---- Motifs (cisBP) and chromVAR deviations --------------------------------------
proj <- addMotifAnnotations(proj, motifSet = "cisbp", name = "Motif")
proj <- addBgdPeaks(proj)
proj <- addDeviationsMatrix(proj, peakAnnotation = "Motif", force = TRUE)   # ADDED: creates MotifMatrix
proj <- saveArchRProject(proj)
