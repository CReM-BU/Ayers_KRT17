# =============================================================================
# config.R — shared paths for all R analysis scripts
#
# Every script in 02_scRNA_analysis/ and 03_scATAC_analysis/ starts with
#   source("config.R")
# and must be run from the repository root.
#
# Paths default to the locations used on the BU SCC. To run elsewhere, set the
# environment variables below (e.g. in ~/.Renviron) instead of editing scripts.
# =============================================================================

# Directory holding the Seurat objects (.rds) and the ArchR project
DATA_DIR <- Sys.getenv(
  "KRT17_DATA_DIR",
  "/restricted/projectnb/crem-trainees/Kotton_Lab/Lauren/25_06_24_lauren_RDS"
)

# Where figure PDFs and tables are written
PLOT_DIR  <- Sys.getenv("KRT17_PLOT_DIR",  file.path(DATA_DIR, "plots"))
TABLE_DIR <- Sys.getenv("KRT17_TABLE_DIR", file.path(DATA_DIR, "tables"))

# Working directory for SCENIC+ inputs/outputs
SCENIC_DIR <- Sys.getenv(
  "KRT17_SCENIC_DIR",
  "/restricted/projectnb/crem-trainees/Kotton_Lab/Lauren/scenic"
)

# ---- Object names (one lineage; see README "Data objects") -------------------
RNA_RAW_RDS       <- file.path(DATA_DIR, "RNA.rds")             # output of 01_scRNA_preprocessing (combined)
RNA_ANNOTATED_RDS <- file.path(DATA_DIR, "rna_annotated.rds")   # output of 02_scRNA_analysis/01_annotate_clusters.R
ARCHR_PROJECT     <- file.path(DATA_DIR, "LJA_250624_ATAC")     # ArchR project directory

for (d in c(PLOT_DIR, TABLE_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- Shared plotting constants ------------------------------------------------
# Cell-type colours. Order matches levels of `subcluster_names` set in
# 01_annotate_clusters.R; colours are named so order no longer matters.
CLUSTER_COLS <- c(
  "iAT1"                     = "brown3",
  "SOX high TGFB-treated"    = "brown1",
  "KRT5-/KRT17+"             = "darkgoldenrod2",
  "Pre KRT5-/KRT17+"         = "darkolivegreen4",
  "Proliferating (distal)"   = "aquamarine3",
  "iAT2"                     = "cyan4",
  "iBC 1"                    = "deepskyblue3",
  "iBC 2"                    = "deepskyblue4",
  "Differentiating Airway"   = "darkorchid4",
  "Proliferating (proximal)" = "deeppink3",
  "PNEC-like"                = "azure3",
  "Non-lung contaminant"     = "azure4"
)

SAMPLE_COLS <- c(
  "iAT1s L+DCI"       = "#D51F26",
  "iAT2s CK+DCI"      = "#272E6A",
  "iBCs PneumacultEx" = "#208A42",
  "iAT2s TGFB+DCI"    = "#89288F"
)

# UMAP theme without axis text/ticks (used for all UMAP panels)
theme_umap <- function() {
  ggplot2::theme(
    axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank()
  )
}

# Seurat-like theme for ArchR embeddings
theme_seurat <- function() {
  ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line  = ggplot2::element_line(colour = "black"),
      axis.text  = ggplot2::element_text(color = "black"),
      legend.key = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank()
    ) +
    theme_umap()
}

plot_path  <- function(name) file.path(PLOT_DIR, name)
table_path <- function(name) file.path(TABLE_DIR, name)
