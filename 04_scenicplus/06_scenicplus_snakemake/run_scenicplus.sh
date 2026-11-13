#!/bin/bash
# =============================================================================
# 04_scenicplus/06_scenicplus_snakemake/run_scenicplus.sh
#
# Runs the full SCENIC+ pipeline (data preparation, motif enrichment with the
# custom cisTarget database, TF/region-to-gene inference, eRegulon building,
# AUCell) with the Snakemake workflow shipped with SCENIC+.
# This replaces the old prepare_data-only command (05_scenicplus_prepare_data.sh).
#
# Inputs (set in config.yaml -- edit the absolute paths before running):
#   ATAC_062425_cistopicobj_30topics_idents.pkl, rna.fixed.h5ad, region_sets/
#     (05_prep_scenicplus_inputs.py)
#   ATAC_062425.regions_vs_motifs.{rankings,scores}.feather
#     (03_create_cistarget_database.sh)
#   motifs-v10-nr.hgnc-m0.00001-o0.0.tbl (aertslab motif annotations)
# Output: scplusmdata.h5mu -> 07_plot_regulon_heatmaps.py
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"

WORKDIR="${WORKDIR:-scplus_pipeline}"

# Create the Snakemake workflow skeleton, then use our config
scenicplus init_snakemake --out_dir "${WORKDIR}"
cp config.yaml "${WORKDIR}/Snakemake/config/config.yaml"

cd "${WORKDIR}/Snakemake"
snakemake --cores 10    # n_cpu: 10 in config.yaml
