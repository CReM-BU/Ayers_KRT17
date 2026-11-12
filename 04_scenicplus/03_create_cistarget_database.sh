#!/bin/bash
# =============================================================================
# 04_scenicplus/03_create_cistarget_database.sh
#
# Step 2 of building the custom cisTarget database: score all aertslab v10
# (non-redundant, clustered) motifs against the padded peak sequences.
#
# Uses create_cistarget_motif_databases.py from
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
MOTIF_LIST="${OUT_DIR}/motifs.txt"   # TODO: document how motifs.txt was made (e.g. ls $CBDIR | sed 's/.cb$//')
DATABASE_PREFIX="ATAC_062425"

"${SCRIPT_DIR}/create_cistarget_motif_databases.py" \
    -f "${FASTA_FILE}" \
    -M "${CBDIR}" \
    -m "${MOTIF_LIST}" \
    -o "${OUT_DIR}/${DATABASE_PREFIX}" \
    --bgpadding 1000 \
    -t 10
