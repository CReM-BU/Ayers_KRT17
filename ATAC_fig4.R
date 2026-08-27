setwd('/restricted/projectnb/crem-trainees/Kotton_Lab/Lauren/25_06_24_lauren_RDS/')
library(ArchR)
library(Seurat, lib.loc = "/share/pkg.8/r/4.4.0/install/lib64/R/library")

# frag_1 <- '/restricted/projectnb/crem-bioinfo/project_workspace/25_06_24_lauren/ATAC/calculations/cellranger_count/AL_KD_ATAC_1/outs/fragments.tsv.gz'
# frag_2 <- '/restricted/projectnb/crem-bioinfo/project_workspace/25_06_24_lauren/ATAC/calculations/cellranger_count/AL_KD_ATAC_2/outs/fragments.tsv.gz'
# frag_3 <- '/restricted/projectnb/crem-bioinfo/project_workspace/25_06_24_lauren/ATAC/calculations/cellranger_count/AL_KD_ATAC_3/outs/fragments.tsv.gz'
# frag_4 <- '/restricted/projectnb/crem-bioinfo/project_workspace/25_06_24_lauren/ATAC/calculations/cellranger_count/AL_KD_ATAC_4/outs/fragments.tsv.gz'
# input_files <- c(frag_1, frag_2, frag_3, frag_4)
# names(input_files) <- c('iAT2', 'iAT1', 'TGFB_iAT2', 'iBC')
# addArchRGenome('hg38')
# addArchRThreads(threads = 1)
# ArrowFiles <- createArrowFiles(inputFiles = input_files, sampleNames = names(input_files)) #min TSS=4 and min frags = 1000
# doubScores <- addDoubletScores(input = ArrowFiles,k = 10,  knnMethod = "UMAP",LSIMethod = 1)
# proj <- ArchRProject(ArrowFiles = ArrowFiles, outputDirectory = "LJA_250624_ATAC", copyArrows = TRUE)
# proj <- filterDoublets(proj)
# proj <- addIterativeLSI(proj, useMatrix = "TileMatrix", name = "IterativeLSI")
# proj <- addClusters(proj, reducedDims = "IterativeLSI")
# saveArchRProject(proj)
# proj <- addUMAP(proj, reducedDims = "IterativeLSI")
# pdf('archr_atac_umap.pdf')
# plotEmbedding(proj)
# plotEmbedding(proj, name = 'Clusters')
# dev.off()
# proj <- loadArchRProject("LJA_250624_ATAC")
# proj <- addGroupCoverages(proj)
# saveArchRProject(proj)
# proj <- addReproduciblePeakSet(proj, groupBy = "Clusters", pathToMacs2 = findMacs2())
# proj <- addPeakMatrix(proj)
# saveArchRProject(proj)
# proj<- addMotifAnnotations(proj, motifSet = "cisbp", name = "Motif") 
# saveArchRProject(proj)




### RNA integration
proj <- loadArchRProject("LJA_250624_ATAC")
rna <- readRDS('rna.subclustered.rds')
addArchRThreads(threads = 2)
DefaultAssay(rna) <- 'RNA'
rna <- JoinLayers(rna)
proj<- addGeneIntegrationMatrix(ArchRProj = proj,  useMatrix = "GeneScoreMatrix",matrixName = "GeneIntegrationMatrix", 
                                 groupRNA = 'subcluster_names',  reducedDims = "IterativeLSI", seRNA = rna, addToArrow = T, force = T)
proj <- saveArchRProject(proj)

plotEmbedding(proj, name = 'predictedGroup') + plotEmbedding(proj, name = 'Sample') + plotEmbedding(proj, name = 'predictedScore') 
proj <- saveArchRProject(proj)


### Clustering of ATAC
proj <- addClusters(
  input = proj,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters_0.25",
  resolution = 0.25
)


### Heatmap to look at label assigments
cM <- confusionMatrix(paste0(proj$Clusters_0.25), paste0(proj$predictedGroup))
library(pheatmap)
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
  mat = as.matrix(cM), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black"
)
## Fig S4B
pdf('plots/ATAC_clusteridents_heatmap.pdf', width = 6, height = 4)
p
dev.off()

## Subclustering of C7 to resolve KRT5-/KRT17+ and Pre KRT5-/KRT17+
idxSample <- BiocGenerics::which(proj$Clusters_0.25 %in% "C7")
projC7 <- subsetArchRProject(
  ArchRProj = proj,
  cells = proj$cellNames[idxSample],
  outputDirectory = "archrC7sub",
  dropCells = TRUE,
  force = TRUE
)

projC7 <- addIterativeLSI(projC7, useMatrix = "TileMatrix", name = "IterativeLSI_C7", force = T)
projC7 <- addClusters(projC7, reducedDims = "IterativeLSI_C7", name = "C7_subclusters", resolution = 0.25, force =T)
saveArchRProject(projC7)

projC7$C7_subclusters <- paste('C7_', projC7$C7_subclusters, sep = '')
projC7$C7_subclusters <- str_replace(projC7$C7_subclusters, "_C", "_")


parent_cells <- getCellNames(proj)
base_vec <- getCellColData(proj, "Clusters_0.25")[,1]
names(base_vec) <- parent_cells

sub_vec <- getCellColData(projC7, "C7_subclusters")[,1]
names(sub_vec) <- getCellNames(projC7)

common <- intersect(parent_cells, names(sub_vec))
merged <- base_vec
merged[common] <- sub_vec[common]

proj <- addCellColData(proj, data = merged, name = "subclusters",
                       cells = names(merged), force = TRUE)
  
##Fig S4A and S4C
pdf('plots/ATAC_clusters.pdf', width = 6, height = 5)
plotEmbedding(proj, name = 'Clusters_0.25', discreteSet = 'calm') + theme_seurat() + xlab('umap_1') + ylab('umap_2') + ggtitle('')
dev.off()
pdf('plots/ATAC_subclusters.pdf', width = 6, height = 5)
plotEmbedding(proj, name = 'subclusters', discreteSet = 'calm') + theme_seurat() + xlab('umap_1') + ylab('umap_2') + ggtitle('')
dev.off()

## Now look at predicted cluster labels with C7 subcluster resolution
cM <- confusionMatrix(paste0(proj$subclusters), paste0(proj$predictedGroup))
library(pheatmap)
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
  mat = as.matrix(cM), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black"
)

## Figure 4B
pdf('plots/ATAC_subclusteridents_heatmap.pdf', width = 6, height = 4)
p
dev.off()




labels <- unique(proj$subclusters)
newLabels <- c("iAT1", 'PNEC-like', 'SOX high TGFB-treated', 'KRT5-/KRT17+', 'iAT2', 'iBC 2', 'iBC 1', 'Pre KRT5-/KRT17+', 'Differentiating Airway')
proj$ClusterLabels <- mapLabels(proj$subclusters,  newLabels = newLabels, oldLabels = labels)

saveArchRProject(proj)
saveArchRProject(projC7)


## Fig 4C: Celltype marker peak heatmap
markerPeaks <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix", 
  groupBy = "ClusterLabels",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon")

markerList <- getMarkers(markerPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")

heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markerPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)

pdf('plots/ATAC-peak-heatmap.pdf', width = 8, height = 6)
heatmapPeaks
dev.off()


## Motif enrichment heatmap (by celltype)
enrichMotifs <- peakAnnoEnrichment(
  seMarker = markerPeaks,
  ArchRProj = proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 20, transpose = TRUE)
pdf('plots/ATAC-motif-heatmap.pdf', width = 8, height = 6)
heatmapEM
dev.off()


## Motif deviations and seqlogos 
proj <- addBgdPeaks(proj)
plotVarDev <- getVarDeviations(proj, name = "MotifMatrix", plot = TRUE)


motifMatches <- getMatches(ArchRProj = proj, name = "Motif")
getPeakSet(proj)[motifMatches@assays@data$matches[,grep('NKX21',colnames(motifMatches))]]

motifs <- c("NKX21", "SMAD3", "TEAD4", "TEAD3", "TEAD1", "TP63", "BACH2", "JUN", 'FOS', 'NFKB1', 'NFKB2')
markerMotifs <- getFeatures(proj, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs <- grep("z:", markerMotifs, value = TRUE)
markerMotifs 
cluster_cols <- c("brown3", "brown1", "darkgoldenrod2", "darkolivegreen4", "aquamarine3", "cyan4", "deepskyblue3", "deepskyblue4", "darkorchid4", "deeppink3", "azure3", "azure4")

names(cluster_cols) <- c( "iAT1", "SOX high TGFB-treated",  "KRT5-/KRT17+", "Pre KRT5-/KRT17+",  "Proliferating (distal)",  "iAT2",  "iBC 1",
  "iBC 2", "Differentiating Airway", "Proliferating (proximal)", "PNEC-like", "Non-lung contaminant")
p <- plotGroups(ArchRProj = proj[proj$ClusterLabels %in% c('Pre KRT5-/KRT17+', 'KRT5-/KRT17+', 'iAT1', 'iAT2', 'iBC 1', 'iBC 2', 'Secretory-like'),], 
                groupBy = "ClusterLabels", 
                colorBy = "MotifMatrix", 
                name = markerMotifs,
                pal = cluster_cols,
                
)
pdf('plots/NKX_deviations.pdf', height = 3, width = 3)
p$`z:NKX21_476`
dev.off()
pdf('plots/SMAD3_deviations.pdf', height = 3, width = 3)
p$`z:SMAD3_743`
dev.off()
pdf('plots/TEAD4_deviations.pdf', height = 3, width = 3)
p$`z:TEAD4_797`
dev.off()
pdf('plots/TP63_deviations.pdf', height = 3, width = 3)
p$`z:TP63_704`
dev.off()
pdf('plots/NFKB1_deviations.pdf', height = 3, width = 3)
p$`z:NFKB1_719`
dev.off()
pdf('plots/FOSL1_deviations.pdf', height = 3, width = 3)
p$`z:FOSL1_142`
dev.off()


theme_seurat <- function() {
  theme_bw(base_size = 14) +
    theme(
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.text = element_text(color = "black"),
      legend.key = element_blank(),
      panel.border = element_blank(),
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(),
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
}
pdf('plots/NKX_motif_umap.pdf', height = 6, width = 5)
plotEmbedding(proj, colorBy = "motifMatrix", name = "z:NKX21_476") + theme_seurat() + xlab('umap_1') + ylab('umap_2') + ggtitle('')
dev.off()
pdf('plots/SMAD3_motif_umap.pdf', height = 6, width = 5)
plotEmbedding(proj, colorBy = "motifMatrix", name = "z:SMAD3_743") + theme_seurat() + xlab('umap_1') + ylab('umap_2') + ggtitle('')
dev.off()
pdf('plots/TEAD4_motif_umap.pdf', height =6, width = 5)
plotEmbedding(proj, colorBy = "motifMatrix", name = "z:TEAD4_797") + theme_seurat() + xlab('umap_1') + ylab('umap_2') + ggtitle('')
dev.off()
pdf('plots/TP63_motif_umap.pdf', height = 6, width = 5)
plotEmbedding(proj, colorBy = "motifMatrix", name = "z:TP63_704") + theme_seurat() + xlab('umap_1') + ylab('umap_2') + ggtitle('')
dev.off()
pdf('plots/NFKB1_motif_umap.pdf', height = 6, width = 5)
plotEmbedding(proj, colorBy = "motifMatrix", name = "z:NFKB1_719") + theme_seurat() + xlab('umap_1') + ylab('umap_2') + ggtitle('')
dev.off()
pdf('plots/FOSL1__motif_umap.pdf', height = 6, width = 5)
plotEmbedding(proj, colorBy = "motifMatrix", name = "z:FOSL1_142") + theme_seurat() + xlab('umap_1') + ylab('umap_2') + ggtitle('')
dev.off()


library(ggseqlogo)
pwm <- getPeakAnnotation(proj, "Motif")$motifs[["FOSL1_142"]]
PWMatrixToProbMatrix <- function(x){
  if (class(x) != "PWMatrix") stop("x must be a TFBSTools::PWMatrix object")
  m <- (exp(as(x, "matrix"))) * TFBSTools::bg(x)/sum(TFBSTools::bg(x))
  m <- t(t(m)/colSums(m))
  m}
ppm <- PWMatrixToProbMatrix(pwm)
p
ggseqlogo(ppm, method = "bits")
