#!/bin/bash
# =============================================================================
# 04_scenicplus/03_create_cistarget_database.sh
#
# Step 2 of building the custom cisTarget database: score all aertslab v10
# (non-redundant, clustered) motifs against the padded peak sequences.
#
# pycistarget 1.1. Uses create_cistarget_motif_databases.py from
#   https://github.com/aertslab/create_cisTarget_databases
# Motif collection: https://resources.aertslab.org/cistarget/motif_collections/
#
# Input : ATAC_062425_padded.fa (from 02_prepare_padded_fasta.sh), motifs.txt
# Output: ATAC_062425.regions_vs_motifs.{rankings,scores}.feather
# =============================================================================
set -euo pipefail

# ---- Edit these paths ----------------------------------------------------------
OUT_DIR="${OUT_DIR:-$PWD}"
SCRIPT_DIR="${SCRIPT_DIR:-/Users/ljayers/Desktop/scenic/create_cisTarget_databases}"
CBDIR="${OUT_DIR}/aertslab_motif_colleciton/v10nr_clust_public/singletons"
FASTA_FILE="${OUT_DIR}/ATAC_062425_padded.fa"
MOTIF_LIST="${OUT_DIR}/motifs.txt"

# motifs.txt = every motif file in the aertslab v10 non-redundant clustered
# collection (singletons), downloaded as described in the pycistarget docs:
#   https://pycistarget.readthedocs.io/en/stable/motif_table.html
# The original zip link is no longer active; the collection is listed at
#   https://motifcollections.aertslab.org/v10nr_clust_public/
# TODO (Lauren): add the motifs.txt used for the paper to this folder.
ls "${CBDIR}" > "${MOTIF_LIST}"
DATABASE_PREFIX="ATAC_062425"

"${SCRIPT_DIR}/create_cistarget_motif_databases.py" \
    -f "${FASTA_FILE}" \
    -M "${CBDIR}" \
    -m "${MOTIF_LIST}" \
    -o "${OUT_DIR}/${DATABASE_PREFIX}" \
    --bgpadding 1000 \
    -t 10
