# scTOP analysis (Fig 3G, 3H, S3H)

Analysis performed by Eitan Vilker

## Core Package Versions
scTOP v2.0.3
SimilarityHelper and TopObject commit x

## Reference Basis Construction
Reference basis used in analysis was constructed from cell types of two datasets. These were:
1. Guo et al LungMAP Human Lung CellRef AT1, AT2, basal, club, goblet and ciliated cells (downsampled
to 500 cells per type)
2. Habermann et al. KRT5-/KRT17+ cells (485 cells)
   
ANOVA feature selection was then applied to select the top 20% of genes for successfully recapitulating original annotations.
