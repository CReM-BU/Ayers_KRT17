# =============================================================================
# 03_scATAC_analysis/02_integrate_cluster_annotate.R
#
# 1. Unconstrained integration with the annotated scRNA-seq object
#    (label transfer -> `predictedGroup`)
# 2. Louvain clustering at resolution 0.25 (`Clusters_0.25`)
# 3. Sub-clustering of cluster C7 to resolve KRT5-/KRT17+ vs Pre KRT5-/KRT17+
#    (`subclusters`)
# 4. Cell-type labels for each (sub)cluster (`ClusterLabels`)
#
# Input : ARCHR_PROJECT, RNA_ANNOTATED_RDS
# Output: ARCHR_PROJECT (updated), C7 sub-project in DATA_DIR/archrC7sub
#
# Run from repo root:  Rscript 03_scATAC_analysis/02_integrate_cluster_annotate.R
# =============================================================================

source("config.R")
suppressPackageStartupMessages({
  library(ArchR)    # v1.0.3
  library(Seurat)   # v5.0.2
  library(stringr)
})
set.seed(1)
addArchRThreads(threads = 2)

proj <- loadArchRProject(ARCHR_PROJECT)

# ---- 1. Integration with scRNA-seq ----------------------------------------------
# NOTE: original script read 'rna.subclustered.rds'; now uses the single annotated
# object produced by 02_scRNA_analysis/01_annotate_clusters.R.
rna <- readRDS(RNA_ANNOTATED_RDS)
DefaultAssay(rna) <- "RNA"
rna <- JoinLayers(rna)

proj <- addGeneIntegrationMatrix(
  ArchRProj = proj, useMatrix = "GeneScoreMatrix", matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI", seRNA = rna, groupRNA = "subcluster_names",
  addToArrow = TRUE, force = TRUE
)
proj <- saveArchRProject(proj)

# ---- 2. Louvain clustering, resolution 0.25 ---------------------------------------
proj <- addClusters(input = proj, reducedDims = "IterativeLSI", method = "Seurat",
                    name = "Clusters_0.25", resolution = 0.25)

# ---- 3. Sub-cluster C7 -------------------------------------------------------------
idxC7 <- BiocGenerics::which(proj$Clusters_0.25 %in% "C7")
projC7 <- subsetArchRProject(ArchRProj = proj, cells = proj$cellNames[idxC7],
                             outputDirectory = file.path(DATA_DIR, "archrC7sub"),
                             dropCells = TRUE, force = TRUE)
projC7 <- addIterativeLSI(projC7, useMatrix = "TileMatrix", name = "IterativeLSI_C7", force = TRUE)
projC7 <- addClusters(projC7, reducedDims = "IterativeLSI_C7", name = "C7_subclusters",
                      resolution = 0.25, force = TRUE)
projC7$C7_subclusters <- str_replace(paste0("C7_", projC7$C7_subclusters), "_C", "_")   # e.g. "C7_1"
saveArchRProject(projC7)

# Merge sub-cluster labels back into the parent project
merged <- setNames(getCellColData(proj, "Clusters_0.25")[, 1], getCellNames(proj))
sub_vec <- setNames(getCellColData(projC7, "C7_subclusters")[, 1], getCellNames(projC7))
common <- intersect(names(merged), names(sub_vec))
merged[common] <- sub_vec[common]
proj <- addCellColData(proj, data = merged, name = "subclusters", cells = names(merged), force = TRUE)

# ---- 4. Cell-type labels -------------------------------------------------------------
# Each (sub)cluster is named after the majority RNA-transferred label (see Fig S4B/D).
# Explicit mapping (same as 04_scenicplus/05_prep_scenicplus_inputs.py and the
# "C#: label" legend of Fig 4B).
label_map <- c(
  C1   = "iBC 2",
  C2   = "iBC 1",
  C3   = "PNEC-like",
  C4   = "iAT1",
  C5   = "SOX high TGFB-treated",
  C6   = "iAT2",
  C7_1 = "KRT5-/KRT17+",
  C7_2 = "Pre KRT5-/KRT17+",
  C8   = "Differentiating Airway"
)
labels <- sort(unique(proj$subclusters))
stopifnot(setequal(labels, names(label_map)))
newLabels <- unname(label_map[labels])
proj$ClusterLabels <- unname(label_map[proj$subclusters])

# Sanity check: assigned label vs. majority predicted RNA label
majority <- tapply(proj$predictedGroup, proj$subclusters, function(x) names(which.max(table(x))))
check <- data.frame(subcluster = labels, assigned = newLabels, majority_RNA_label = majority[labels])
print(check)
if (any(check$assigned != check$majority_RNA_label)) {
  warning("Some assigned labels differ from the majority RNA label -- check the table above.")
}

saveArchRProject(proj)
