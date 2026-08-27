library(Seurat)
library(ggplot2)
library(presto)
library(dplyr)
library(RColorBrewer)

setwd('/restricted/projectnb/crem-trainees/Kotton_Lab/Lauren/25_06_24_lauren_RDS')
#atac <- readRDS('ATAC.rds')
rna <- readRDS('RNA.rds')

DimPlot(rna, group.by = "sample")+
DimPlot(rna, group.by = "SCT_snn_res.0.25", label = T)
DimPlot(rna, group.by = "sample")+
DimPlot(rna, group.by = "SCT_snn_res.0.5", label = T)+
  FeaturePlot(rna, c('MKI67', 'TOP2A'))


meta.data <- rna[[]]
counts <- group_by(meta.data, sample, SCT_snn_res.0.25) %>% summarise(count = n())
ggplot(counts, aes(SCT_snn_res.0.25, count, fill = sample)) +
  geom_bar(stat = 'identity', position = 'fill')

DefaultAssay(rna) <- 'SCT'
rna$cluster_names <- ifelse(rna$SCT_snn_res.0.5 %in% c(10, 0, 9), 'iAT1',
                            ifelse(rna$SCT_snn_res.0.5 %in% c(8,6), 'iAT2',
                                   ifelse(rna$SCT_snn_res.0.5 == 5, 'TGFB 1',
                                          ifelse(rna$SCT_snn_res.0.5 == 2, 'TGFB2',
                                                 ifelse(rna$SCT_snn_res.0.5 == 14, 'TGFB or LDCI',
                                                        ifelse(rna$SCT_snn_res.0.5 %in% c(7,3), 'proliferating 1',
                                                               ifelse(rna$SCT_snn_res.0.5 == 11, 'proliferating 2',
                                                                      paste('iBC', rna$SCT_snn_res.0.5, sep = '-'))))))))

rna <- JoinLayers(rna, assay = 'RNA')
rna <- PrepSCTFindMarkers(rna)
Idents(rna) <- "cluster_names"
sample.markers <- FindAllMarkers(rna, only.pos = TRUE)
sample.markers %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1)
# saveRDS(rna, 'RNA-updated.rds')
rna <- readRDS('RNA-updated.rds')
top20 <- sample.markers %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1) %>% slice_head(n = 20) %>%  ungroup() 
top50 <- sample.markers %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1) %>% slice_head(n = 50) %>%  ungroup() 
write.table(top50, 'cluster_names_top50_markers.tsv', sep = '\t', quote =F)
write.table(sample.markers, 'cluster_names_allmarkers.tsv', sep = '\t', quote =F)
sample.markers <- read.table('cluster_names_allmarkers.tsv', sep = '\t')
pdf('sample_cluster_names_UMAP.pdf', width =14, height =5)
DimPlot(rna, group.by = "sample")+
  DimPlot(rna, group.by = "SCT_snn_res.0.5", label = T)+
  DimPlot(rna, group.by = "cluster_names", label = T)
dev.off()


Idents(rna) <- "SCT_snn_res.1"
new.cluster.ids <- c("iBC", "iAT1", "iBC", "unk_TGFB_1", "iAT2", "ABI 2", "iAT1", "distal proliferating", "iAT2", "iBC", 
                     "proximal proliferating", "distal proliferating", "unk_iAT1", "distal proliferating", "iAT1", "ABI 1",
                     "iBC", "unk_TGFB_2", "iBC", "TGFB_3", "unk_iBC")
names(new.cluster.ids) <- levels(rna)
rna <- RenameIdents(rna, new.cluster.ids)
DimPlot(rna, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
rna$cluster_names <- Idents(rna)
#saveRDS(rna, "RNA_082825.rds")

rna <- readRDS("RNA_082825.rds")

AT1.AT2 <- FindMarkers(rna, ident.1 = "iAT1", ident.2 = "iAT2")
