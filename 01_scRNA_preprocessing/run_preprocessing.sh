#!/bin/bash
# =============================================================================
# 01_scRNA_preprocessing/run_preprocessing.sh
#
# Renders the per-sample and combined Seurat pipelines.
#
# Reconstructed from the Rmd params blocks (the render commands themselves were
# not saved). Mitochondrial cut-off: cells with >= 25% mito reads are removed
# in 10x.Rmd (per sample). The combined Rmd does not filter on percent.mito.
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"

# ---- Per-sample QC / processing --------------------------------------------------
PERCENT_MITO=25
for pref in KD_AL_OCM_1 KD_AL_OCM_2 KD_AL_OCM_3 KD_AL_OCM_4; do
  Rscript -e "rmarkdown::render('10x.Rmd', output_file = '${pref}.html',
    params = list(prefix = '${pref}', resDEG = 'SCT_snn_res.0.25',
                  percent.mito = ${PERCENT_MITO}, cell.source = 'human'))"
done

# ---- Combined object ---------------------------------------------------------------
Rscript -e "rmarkdown::render('10x.combined.Rmd', output_file = 'combined.html',
  params = list(prefix = 'combined', resDEG = 'SCT_snn_res.0.25', percent.mito = 0.25,
                cell.source = 'human', harmonize = FALSE,
                sample.names = c('KD_AL_OCM_1','KD_AL_OCM_2','KD_AL_OCM_3','KD_AL_OCM_4')))"
