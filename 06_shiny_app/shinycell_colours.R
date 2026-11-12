# =============================================================================
# 06_shiny_app/shinycell_colours.R
#
# Colour / label settings for the ShinyCell app
#   https://crem-bu.shinyapps.io/scRNAseq-TGFB-treatment-iAT2s/
# from which the published Fig 3D (GFP, KRT17) and Fig S3E/S3F (ASCL1, CHGA,
# SOX2/4/6) feature plots were exported.
#
# Source: 25_06_25_lauren_RNA-cluster_annotations.R (L. Ayers)
#
# PARTIAL: the ShinyCell build itself (createConfig() / makeShinyApp()) is not
# yet in the repo. `seu` is the annotated object (RNA_ANNOTATED_RDS) and
# `scConf` comes from ShinyCell::createConfig(seu).
# Note the app uses older cluster names (ABI1 / ABI2 / NKX2.1 negative).
# =============================================================================

library(ShinyCell)

# App labels (paper name -> name used in the app)
seu$Cluster_annotation <- plyr::mapvalues(
  seu$subcluster_names,
  c("Differentiating Airway", "iAT1", "iAT2", "iBC 1", "iBC 2", "KRT5-/KRT17+",
    "Non-lung contaminant", "PNEC-like", "Pre KRT5-/KRT17+", "Proliferating (distal)",
    "Proliferating (proximal)", "SOX high TGFB-treated"),
  c("Differentiating Airway", "iAT1", "iAT2", "iBC 1", "iBC 2", "ABI2",
    "NKX2.1 negative", "PNEC-like", "ABI1", "Proliferating (distal)",
    "Proliferating (proximal)", "SOX high TGFB-treated")
)
# NOTE: the original line mapped from seu$TGFB.subclustered (cluster IDs such as
# "1_1"), which would not match these names; subcluster_names is used here.

scConf <- modColours(scConf, meta.to.mod = "Cluster_annotation",
                     new.colours = c("darkolivegreen4", "darkgoldenrod2", "darkorchid4", "brown3",
                                     "cyan4", "deepskyblue3", "deepskyblue4", "azure4", "azure3",
                                     "aquamarine3", "deeppink3", "brown1"))
scConf <- modColours(scConf, meta.to.mod = "orig.ident",
                     new.colours = c("#D51F26", "#272E6A", "#89288F", "#208A42"))
