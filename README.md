# Ayers_KRT17

Code for the single-cell analyses in:

> Ayers L.J.\*, Dooley P.E.\*, *et al.* **Human iPSC modeling of pulmonary alveolar cell fates reveals the aberrant KRT5-/KRT17+ state is profibrotic yet retains alveolar type 1 cell competence.** (manuscript in preparation)

This repository contains the scRNA-seq, scATAC-seq, SCENIC+ and scTOP analyses
behind **Figures 3, 4, S3 and S4**. See **[FIGURES.md](FIGURES.md)** for the
panel-by-panel map of script → output file.

Interactive scRNA-seq browser: https://crem-bu.shinyapps.io/scRNAseq-TGFB-treatment-iAT2s/

## Data

Samples (BU1 SFTPC-tdTomato / KRT17-GFP iPSC line, clone 42):

| Sample | scRNA-seq (`sample`) | scATAC-seq (ArchR `Sample`) |
|---|---|---|
| iAT2s in CK+DCI | `iAT2s CK+DCI` | `iAT2` |
| iAT2s, 48 h TGF-β in DCI | `iAT2s TGFB+DCI` | `TGFB_iAT2` |
| iAT1s in L+DCI | `iAT1s L+DCI` | `iAT1` |
| iBCs in Pneumacult Ex+ | `iBCs PneumacultEx` | `iBC` |

Raw and processed data: GEO accession *TBD*.

## Repository layout

Run the steps in order. All R scripts are run from the repository root and read
paths from [`config.R`](config.R).

```
00_alignment/              FASTQ generation and alignment (BU CReM bioinformatics core, SCC)
  scRNA/                   reference build, cellranger 9.0.1 mkfastq + multi (OCM)
  scATAC/                  cellranger-atac 2.1.0 mkfastq + count
01_scRNA_preprocessing/    per-sample QC + merged Seurat object (Rmd)
02_scRNA_analysis/
  01_annotate_clusters.R   sub-cluster cluster 1, name cell types -> rna_annotated.rds
  02_markers_and_DEG.R     MAST cluster markers, KRT5-/KRT17+ vs iAT2
  fig3.R, figS3.R          Figure 3 and S3 panels
03_scATAC_analysis/
  01_build_archr_project.R QC, doublets, LSI, peaks, motifs, chromVAR
  02_integrate_cluster_annotate.R  RNA label transfer, clustering, C7 sub-clusters
  fig4.R, figS4.R          Figure 4 and S4 panels
04_scenicplus/             eRegulon analysis (pycisTopic + custom cisTarget db + SCENIC+)
05_scTOP/                  scTOP projections (pending)
06_shiny_app/              ShinyCell app build (Fig 3D, S3E, S3F were exported from it)
config.R                   paths, colours, themes shared by all R scripts
```

## Data objects

One object lineage is used throughout (paths set in `config.R`):

```
cellranger multi ──> 01_scRNA_preprocessing ──> combined/rds/sc.Rds (shared as RNA.rds)
                                                   │ 02_scRNA_analysis/01_annotate_clusters.R
                                                   ▼
                                           rna_annotated.rds ──> fig3.R, figS3.R
                                                   │                    ▲
cellranger-atac ──> 03/01_build_archr_project.R    │            02_markers_and_DEG.R
                         │                         ▼
                         └──> 03/02_integrate_cluster_annotate.R ──> fig4.R, figS4.R
                                                   │
                                                   ▼
                                       04_scenicplus/01_export_inputs.R ──> SCENIC+
```

## Software

**scRNA-seq preprocessing** (`01_scRNA_preprocessing`; from `sessionInfo()`, BU SCC, AlmaLinux 8.10)

| Tool | Version |
|---|---|
| R | 4.4.3 |
| Seurat / SeuratObject | 5.5.1 / 5.4.0 |
| sctransform | 0.4.3 |
| glmGamPoi | 1.18.0 |
| presto | 1.0.0 |
| tidyverse / dplyr / data.table | 2.0.0 / 1.2.1 / 1.17.8 |
| clustree | 0.5.1 |

**Alignment**

| Tool | Version |
|---|---|
| cellranger mkref (scRNA reference) | 3.0.2 |
| cellranger mkfastq / multi (scRNA) | 9.0.1 |
| cellranger-atac | 2.1.0 |

**Downstream analysis and figures** (`02`–`06`)

| Tool | Version |
|---|---|
| Seurat | 5.0.2 *(lib.loc in original figure scripts; to confirm)* |
| ArchR | 1.0.3 |
| MAST, fgsea, msigdbr, pheatmap, ggseqlogo | *TBD*: printed by `sessionInfo()` at the end of each figure script |
| pycisTopic, SCENIC+ | *TBD* |
| scTOP | 2.0.3 |
| ShinyCell | *TBD* |

## Status

Code still to be added is listed in [FIGURES.md](FIGURES.md) (❌ / ⏳) and in the
`TODO` comments inside scripts.
