#!/bin/bash
# =============================================================================
# 01_scRNA_preprocessing/run_preprocessing.sh
#
# Renders the per-sample and combined Seurat pipelines.
#
# RECONSTRUCTED from the Rmd params blocks -- the exact render commands used
# for the paper were not in the repo. CONFIRM percent.mito for each sample
# (methods state 15-35%; the Rmd default is 25).
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"

# ---- Per-sample QC / processing --------------------------------------------------
declare -A PERCENT_MITO=(          # TODO: confirm values used for each sample
  [KD_AL_OCM_1]=25
  [KD_AL_OCM_2]=25
  [KD_AL_OCM_3]=25
  [KD_AL_OCM_4]=25
)
for pref in KD_AL_OCM_1 KD_AL_OCM_2 KD_AL_OCM_3 KD_AL_OCM_4; do
  Rscript -e "rmarkdown::render('10x.Rmd', output_file = '${pref}.html',
    params = list(prefix = '${pref}', resDEG = 'SCT_snn_res.0.25',
                  percent.mito = ${PERCENT_MITO[$pref]}, cell.source = 'human'))"
done

# ---- Combined object ---------------------------------------------------------------
Rscript -e "rmarkdown::render('10x.combined.Rmd', output_file = 'combined.html',
  params = list(prefix = 'combined', resDEG = 'SCT_snn_res.0.25', percent.mito = 0.25,
                cell.source = 'human', harmonize = FALSE,
                sample.names = c('KD_AL_OCM_1','KD_AL_OCM_2','KD_AL_OCM_3','KD_AL_OCM_4')))"
