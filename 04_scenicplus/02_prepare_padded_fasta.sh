#!/bin/bash
# =============================================================================
# 04_scenicplus/02_prepare_padded_fasta.sh
#
# Step 1 of building a custom cisTarget database from the ArchR consensus peak
# set: extract peak sequences from hg38 with 1 kb background padding.
#
# Uses create_fasta_with_padded_bg_from_bed.sh from
#   https://github.com/aertslab/create_cisTarget_databases
#
# Input : ArchR_allPeaks.bed (from 01_export_inputs.R)
# Output: ATAC_062425_padded.fa
# =============================================================================
set -euo pipefail

# ---- Edit these paths ----------------------------------------------------------
SCRIPT_DIR="${SCRIPT_DIR:-/Users/ljayers/Desktop/scenic/create_cisTarget_databases}"
REGION_BED="${REGION_BED:-/Users/ljayers/Desktop/scenicoutput/ATAC_06_24_25/region_sets/ArchR_allPeaks.bed}"
GENOME_FASTA="${GENOME_FASTA:-/Users/ljayers/Desktop/scenic/reference/hg38.fa}"
CHROMSIZES="${CHROMSIZES:-/Users/ljayers/Desktop/scenic/reference/hg38.chrom.sizes}"
DATABASE_PREFIX="ATAC_062425"

"${SCRIPT_DIR}/create_fasta_with_padded_bg_from_bed.sh" \
    "${GENOME_FASTA}" \
    "${CHROMSIZES}" \
    "${REGION_BED}" \
    "${DATABASE_PREFIX}_padded.fa" \
    1000 \
    yes
