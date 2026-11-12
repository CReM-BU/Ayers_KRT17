#!/bin/bash
# =============================================================================
# 04_scenicplus/05_scenicplus_prepare_data.sh
#
# SCENIC+ step: combine the cisTopic object (ATAC) and the RNA AnnData into a
# MuData object. Data are not multiome, so cells are paired into metacells by
# cell type (GEX:ident) with --is_not_multiome.
#
# Inputs : *_cistopicobj_30topics_idents.pkl  (pycisTopic object with cell-type
#                                              labels -- step that creates it is MISSING)
#          rna.fixed.h5ad                      (from 01_export_inputs.R -- the "fix"
#                                              step is MISSING)
# Output : ACC_GEX.h5mu
#
# The remaining SCENIC+ steps (region-to-gene, TF-to-gene, motif enrichment,
# eRegulon building) and the eRegulon plotting code for Fig 4D / S4E are MISSING.
# =============================================================================
set -euo pipefail

SCENIC_OUT="${SCENIC_OUT:-/Users/ljayers/Desktop/scenicoutput/ATAC_06_24_25}"
SCENIC_RNA="${SCENIC_RNA:-/Users/ljayers/Desktop/scenic/rna}"
SCENIC_RUN="${SCENIC_RUN:-/Users/ljayers/Desktop/scenic/outs}"

scenicplus prepare_data prepare_GEX_ACC \
    --cisTopic_obj_fname "${SCENIC_OUT}/ATAC_062425_cistopicobj_30topics_idents.pkl" \
    --GEX_anndata_fname  "${SCENIC_RNA}/rna.fixed.h5ad" \
    --out_file           "${SCENIC_RUN}/ACC_GEX.h5mu" \
    --is_not_multiome \
    --key_to_group_by GEX:ident \
    --nr_cells_per_metacells 10
