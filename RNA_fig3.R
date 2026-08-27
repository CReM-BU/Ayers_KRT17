library(Seurat, lib.loc = "/share/pkg.8/r/4.4.0/install/lib64/R/library") #v5.0.2
library(ggplot2)
library(presto)
library(dplyr)
library(RColorBrewer)

setwd('/restricted/projectnb/crem-trainees/Kotton_Lab/Lauren/25_06_24_lauren_RDS')
rna <- readRDS('RNA.rds')


## Fig 3B (UMAP by sample)
col <- c("#D51F26", "#272E6A", "#208A42", "#89288F")
names(col) <- c('iAT1s L+DCI', 'iAT2s CK+DCI', 'iBCs PneumacultEx', 'iAT2s TGFB+DCI')
pdf('plots/rna_UMAP_sample.pdf', width = 10, height = 11)
DimPlot(rna, group.by = 'sample', cols=col, pt.size = 1.5) &   
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(),  legend.position = "none") 
dev.off()


DefaultAssay(rna) <- 'SCT'

## make subclusters
Idents(rna) <- "SCT_snn_res.0.25"
rna <- FindSubCluster(rna, '1', subcluster.name = "TGFB.subclustered", graph.name = 'SCT_nn', resolution = 0.25, algorithm = 1) #only 2 cells in sublcuster 1_3, merge with 1_1
DimPlot(rna, group.by = 'TGFB.subclustered', label = T) 
rna$TGFB.subclustered <- ifelse(rna$TGFB.subclustered == "1_3", "1_1", rna$TGFB.subclustered)
cluster_cols <- c("brown3", "brown1", "darkgoldenrod2", "darkolivegreen4", "aquamarine3", "cyan4", "deepskyblue3", "deepskyblue4", "darkorchid4", "deeppink3", "azure3", "azure4")
Idents(rna) <- "TGFB.subclustered"
names(cluster_cols) <- levels(rna)

#pdf('plots/UMAP.tgfb.subclustered.pdf', width = 10, height = 11)
#DimPlot(rna, group.by = 'TGFB.subclustered', label = T, cols = cluster_cols) 
#dev.off()
Idents(rna)
rna <- RenameIdents(rna, "0" = 'iAT1', "1_0" = "SOX high TGFB-treated", "1_1" = "KRT5-/KRT17+", "1_2" = "Pre KRT5-/KRT17+", "2" = "Proliferating (distal)",
                    "3" = "iAT2", "4" = 'iBC 1', "5" = "iBC 2", '6' = "Differentiating Airway", "7" = "Proliferating (proximal)", "8" = "PNEC-like", "9" = 'NKX2.1 negative')
rna$subcluster_names <- Idents(rna)
names(cluster_cols) <- levels(rna)

saveRDS(rna, 'rna.subclustered.updated.rds')
rna <- readRDS('rna.subclustered.updated.rds')
## Fig 3C (UMAP by cluster)
pdf('plots/UMAP.tgfb.subclustered.cellnames.updated.pdf', width = 10, height = 11)
DimPlot(rna, group.by = 'subcluster_names', label = T, cols = cluster_cols, repel = T, label.size = 8, pt.size = 1.3) +  
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), 
        legend.position = "none")
dev.off()

## Figure 3D (GFP and KRT17 expression) made with Shiny app: https://crem-bu.shinyapps.io/scRNAseq-TGFB-treatment-iAT2s/

## Figure 3E (dotplot)
pdf('plots/cluster.dotplot.updated.022326.pdf', height = 6, width = 17)
DotPlot(rna, features = c('NKX2-1', 
                          'AGER', 'CAV1', 'EMP2',
                          'SFTPC', 'SFTPB', 'ETV5',
                          'ICAM1', 'PTGS2', 'TM4SF1','SFN', 'CLDN4', 'TPM1',
                          'ITGB6', 'COL1A1', 'TAGLN', 'FN1', 'ITGAV', 'SMAD3', 'SOX4', 'TGFB1',
                          'PLAU', 'CDKN2B', 
                          'TP53', 'GADD45A', 'CASP3', 'TNF',
                          'KRT17', 'KRT7', 'KRT8', 'KRT19', 
                          'ITGB4', 'S100A2', 'TP63', 'BCAM', 'KRT5',
                          'SCGB1A1', 'CCNO', 'FOXJ1',
                          'TOP2A', 'MKI67', 'CCNB1'),
  col.min = -1, col.max = 5, dot.scale = 10) + 
  scale_y_discrete(limits = c('iAT1', 'iAT2', 'Pre KRT5-/KRT17+', 'KRT5-/KRT17+', 'SOX high TGFB-treated', 'PNEC-like', 'iBC 1', 'iBC 2', 'Differentiating Airway', 'Proliferating (proximal)', 'Proliferating (distal)')) + 
  RotatedAxis() + xlab('') + ylab('')
dev.off()

## Figure 3F (module violin plots)
pdf('plots/haberman_vlnplots.pdf', width = 12, height = 9)
VlnPlot(rna, features = c('Hb_AT2', 'Hb_KRT5neg.KRT17pos', 'Hb_AT1', 'Hb_Basal', 'Hb_MUC5Bpos'), pt.size = 0, cols = cluster_cols)
dev.off()

lvl <-  c("iAT2", "Pre KRT5-/KRT17+", "KRT5-/KRT17+", "iAT1",
          "iBC 1", "iBC 2", "PNEC-like", "Proliferating (distal)", "SOX high TGFB-treated",
"Proliferating (proximal)", "Differentiating Airway",   "NKX2.1 negative")
Idents(rna) <- factor(x = rna$subcluster_names, levels = lvl)
pdf('plots/haberman_vlnplots_subset.pdf', width = 8, height = 5)
VlnPlot(rna, features = c('Hb_AT2', 'Hb_KRT5neg.KRT17pos'), 
        idents = c('iAT2', 'Pre KRT5-/KRT17+', 'KRT5-/KRT17+', 'iAT1', 'iBC 1', 'iBC 2'), 
        pt.size = 0, cols = cluster_cols)
dev.off()

## Figure 3K (gene expression violin plots)
pdf('plots/geneexp_vlnplots.pdf', width = 12, height = 9)
VlnPlot(rna, features = c('ICAM1', 'TAGLN', 'CDKN1A', 'COL1A1', 'IL32'), pt.size = 0, cols = cluster_cols)
dev.off()
pdf('plots/geneexp_IL32.pdf', width = 6, height = 5)
VlnPlot(rna, features = c('IL32'), pt.size = 0, cols = cluster_cols)
dev.off()

## Figure 3L from shiny app: https://crem-bu.shinyapps.io/19_09_10_darrell_habermann/

### DEGs and Hallmark (Figure 3I)
library(msigdbr)
library(fgsea)

# rna <- PrepSCTFindMarkers(rna)
# Cluster markers (Table XX)
# cluster.markers <- FindAllMarkers(rna, test.use = "MAST", only.pos = F, min.pct = 0.1, verbose = F)
# top20 <- cluster.markers %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1) %>% slice_head(n = 20) %>%  ungroup() 
# top50 <- cluster.markers %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1) %>% slice_head(n = 50) %>%  ungroup() 
# write.table(top20, 'plots/cluster_names_top20_markers.tsv', sep = '\t', quote =F)
# write.table(top50, 'plots/cluster_names_top50_markers.tsv', sep = '\t', quote =F)
# write.table(cluster.markers, 'plots/cluster_names_allmarkers.tsv', sep = '\t', quote =F)


## KRT5-/KRT17+ vs iAT2 (Hallmark)
# ABI2.iAT2.markers <- FindMarkers(rna, ident.1 = 'KRT5-/KRT17+', ident.2 = 'iAT2', test.use = "MAST", only.pos = F, min.pct = 0.1, verbose = F)
# write.table(ABI2.iAT2.markers, 'plots/ABI2.iAT2.markers.sct.tsv', sep = '\t', quote = F, col.names = T, row.names = T)
ABI2.iAT2.markers <- read.table('plots/ABI2.iAT2.markers.sct.tsv', sep = '\t')
ranks <- c(ABI2.iAT2.markers$avg_log2FC)
names(ranks) <- rownames(ABI2.iAT2.markers)
msig <- msigdbr(species = 'Homo sapiens', category = 'H') %>% dplyr::select(gs_name, gene_symbol)
pathways <- split(msig$gene_symbol, msig$gs_name)

fgsea_res <- fgsea(pathways = pathways, stats = ranks, minSize = 5, maxSize = 500)


fgsea_res <- fgsea_res %>%
  arrange(dplyr::desc(abs(NES))) %>%
  mutate(
    Significance = ifelse(padj < 0.05, "padj < 0.05", "Not sig."),
    Significance = factor(Significance, levels = c("padj < 0.05", "Not sig."))
  )
fwrite(fgsea_res, file = "plots/ABI2.iAT2.hallmark.sct.tsv", sep = "\t")

## bar plot
p <- ggplot(fgsea_res, aes(x = reorder(pathway, NES), y = NES, fill = Significance)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal(base_size = 11) +
  ggtitle('Enriched Hallmark Pathways:/n ABI2 vs. iAT2') +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 8, face = "bold"), plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(name = "Significance",values = setNames(c("Dark Blue", "Light Blue"), c("padj < 0.05", "Not sig.") )) +
  ylab("Normalized Enrichment Score (NES)")


library(stringr)
fgsea_res <- fgsea_res %>% 
  transmute(pathway = str_replace_all(pathway, "^HALLMARK_", "") |> str_replace_all("_"," "),
            NES = NES, padj = padj) %>%
  arrange(desc(NES))

top_n <- 50
fgsea_res <- fgsea_res %>% 
  slice_max(order_by = abs(NES), n = top_n) %>%
  mutate(
    pathway = factor(pathway, levels = rev(pathway)),
    mlog10p = -log10(padj),
    sig = padj < 0.05
  )
fgsea_res <- fgsea_res[fgsea_res$sig == T,]
fgsea_res$pathway <- with(fgsea_res, reorder(pathway, NES))

pdf('plots/hallmark_ABIvsiAT2_SCT_ordered.pdf', width = 8, height = 7)
ggplot(fgsea_res, aes(x = NES, y = pathway)) +
  geom_segment(aes(x = 0, xend = NES, yend = pathway), linewidth = 0.6, color = "grey80") +
  geom_point(aes(size = mlog10p, color = NES), shape = 16) +
  scale_color_gradient2(low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0,
                        limits = c(-3, 3), oob = scales::squish, name = "NES") +
  scale_size_area(max_size = 10, name = "-log10(padj)") +
  geom_vline(xintercept = 0, linetype = 3, color = "grey70") +
  labs(x = "Normalized Enrichment Score (NES)",
       y = NULL,
       title = "Enriched HALLMARK Pathways: \n KRT5-/KRT17+ vs. iAT2") +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(size = 12),
        legend.position = "right")
dev.off()


### Fig S3A
pdf('plots/rna_UMAP_0.25clusters.pdf', width = 10, height = 11)

DimPlot(rna, group.by = 'SCT_snn_res.0.25', label = T, repel = T, label.size = 8, pt.size = 1.3) +  
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), 
        legend.position = "none")

dev.off()

### Fig S3B

pdf('plots/rna_UMAP_0.25clusters_withsubclustering.pdf', width = 10, height = 11)

DimPlot(rna, group.by = 'TGFB.subclustered', label = T, repel = T, label.size = 8, pt.size = 1.3) +  
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), 
        legend.position = "none") + ggtitle('SCT_snn_res.0.25 with subclustering')

dev.off()


### Fig S3C and D
# switch to RNA normalized for heatmap visualization
markers <- read.table('cluster_names_allmarkers.tsv', sep = '\t')
DefaultAssay(rna) <- "RNA"
rna <- NormalizeData(rna)
all.genes <- rownames(rna)
rna <- ScaleData(rna, features = all.genes, vars.to.regress="percent.mito", verbose = F) # default: only scales var.genes, which misses some in heatmap
g1 <- DoHeatmap(rna, features = top20[1:120,]$gene) + scale_fill_gradientn(colors = c("blue", "white", "red"))
g2 <- DoHeatmap(rna, features = top20[121:240,]$gene) + scale_fill_gradientn(colors = c("blue", "white", "red"))
pdf('plots/genes_heatmap.pdf', width = 10, height = 15)
g1
g2
dev.off()
top10 <- markers %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1) %>% slice_head(n = 10) %>%  ungroup() 
g <- DoHeatmap(rna, features = top10$gene) + scale_fill_gradientn(colors = c("blue", "white", "red"))

DefaultAssay(sc) <- "SCT"

## Fig S3E
pdf('plots/nkx_vlnplots.pdf', width = 12, height = 9)
VlnPlot(rna, features = c('NKX2-1'), pt.size = 0, cols = cluster_cols)
dev.off()



### Heatmap (Figure 3J)

library(Matrix)
library(pheatmap)

#get KRT5-/KRT17+ top 20
top20 <- read.table('plots/cluster_names_top20_markers.tsv', sep = '\t')
top20 <- top20[top20$cluster == 'KRT5-/KRT17+',]

#get Hb K5-/K17 top 20
k17 <- read.csv('modules/habermann.msig.SPRING.csv')
k17 <- k17[9,c(2:21)]

genes      <- as.character(k17)
cap_z      <- 2              
DefaultAssay(rna) <- 'SCT'
rna_sub <- rna #subset(rna, subset = subcluster_names %in% c('iAT1', 'iAT2', 'Pre KRT5-/KRT17+', 'KRT5-/KRT17+', 'iBC 1', 'iBC 2','SOX high TGFB-treated'))

genes_use <- intersect(genes, rownames(rna_sub))
clusters <- as.character(rna_sub@meta.data[["subcluster_names"]])
names(clusters) <- colnames(rna_sub)
mat <- GetAssayData(rna_sub, assay = "SCT", layer = 'data')[genes_use, , drop = FALSE]

avg <- vapply(sort(unique(clusters)), function(cl) {
  cells <- names(clusters)[clusters == cl]
  if (length(cells) == 1) {
    as.numeric(mat[, cells, drop = FALSE])
  } else {
    Matrix::rowMeans(mat[, cells, drop = FALSE])
  }
}, FUN.VALUE = numeric(length(genes_use)))

avg <- as.matrix(avg)
rownames(avg) <- genes_use
colnames(avg) <- sort(unique(clusters))

#scale rows
mu <- rowMeans(avg, na.rm = TRUE)
sdv <- apply(avg, 1, sd, na.rm = TRUE)
sdv[sdv == 0] <- 1
avg_plot <- (avg - mu) / sdv
avg_plot[avg_plot >  cap_z] <-  cap_z
avg_plot[avg_plot < -cap_z] <- -cap_z


ord <- c('iAT2', 'Pre KRT5-/KRT17+', 'KRT5-/KRT17+', 'iAT1', 'iBC 1', 'iBC 2', 'Differentiating Airway', 'SOX high TGFB-treated', 'PNEC-like', 'Proliferating (distal)', 'Proliferating (proximal)', "Non-lung contaminant")
avg_plot <- avg_plot[,ord]
pheatmap(
  avg_plot,
  cluster_rows = TRUE,
  cluster_cols = F,
  angle_col = 45,
  border_color = NA,
  filename = 'plots/haberman_top20_k5negk17pos_heatmap.pdf',
  height = 5,
  width = 5
)



genes      <- top20$gene
cap_z      <- 2              
DefaultAssay(rna) <- 'SCT'
rna_sub <- rna #subset(rna, subset = subcluster_names %in% c('iAT1', 'iAT2', 'Pre KRT5-/KRT17+', 'KRT5-/KRT17+', 'iBC 1', 'iBC 2','SOX high TGFB-treated'))

genes_use <- intersect(genes, rownames(rna_sub))
clusters <- as.character(rna_sub@meta.data[["subcluster_names"]])
names(clusters) <- colnames(rna_sub)
mat <- GetAssayData(rna_sub, assay = "SCT", layer = 'data')[genes_use, , drop = FALSE]

avg <- vapply(sort(unique(clusters)), function(cl) {
  cells <- names(clusters)[clusters == cl]
  if (length(cells) == 1) {
    as.numeric(mat[, cells, drop = FALSE])
  } else {
    Matrix::rowMeans(mat[, cells, drop = FALSE])
  }
}, FUN.VALUE = numeric(length(genes_use)))

avg <- as.matrix(avg)
rownames(avg) <- genes_use
colnames(avg) <- sort(unique(clusters))

#scale rows
mu <- rowMeans(avg, na.rm = TRUE)
sdv <- apply(avg, 1, sd, na.rm = TRUE)
sdv[sdv == 0] <- 1
avg_plot <- (avg - mu) / sdv
avg_plot[avg_plot >  cap_z] <-  cap_z
avg_plot[avg_plot < -cap_z] <- -cap_z


ord <- c('iAT2', 'Pre KRT5-/KRT17+', 'KRT5-/KRT17+', 'iAT1', 'iBC 1', 'iBC 2', 'Differentiating Airway', 'SOX high TGFB-treated', 'PNEC-like', 'Proliferating (distal)', 'Proliferating (proximal)', "Non-lung contaminant")
avg_plot <- avg_plot[,ord]
pheatmap(
  avg_plot,
  cluster_rows = TRUE,
  cluster_cols = F,
  angle_col = 45,
  border_color = NA,
  filename = 'plots/abi2_top20_heatmap.pdf',
  height = 5,
  width = 7
)





