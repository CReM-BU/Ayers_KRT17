OUT_DIR=""${PWD}""
CBDIR="${OUT_DIR}/aertslab_motif_colleciton/v10nr_clust_public/singletons"
FASTA_FILE="${OUT_DIR}/ATAC_062425_padded.fa"
MOTIF_LIST="${OUT_DIR}/motifs.txt"
SCRIPT_DIR="/Users/ljayers/Desktop/scenic/create_cisTarget_databases"

REGION_BED="/Users/ljayers/Desktop/scenicoutput/ATAC_06_24_25/region_sets/ArchR_allPeaks.bed"
GENOME_FASTA="/Users/ljayers/Desktop/scenic/reference/hg38.fa"
CHROMSIZES="/Users/ljayers/Desktop/scenic/reference/hg38.chrom.sizes"
DATABASE_PREFIX="ATAC_062425"

"${SCRIPT_DIR}/create_cistarget_motif_databases.py" \
    -f ${FASTA_FILE} \
    -M ${CBDIR} \
    -m ${MOTIF_LIST} \
    -o ${OUT_DIR}/${DATABASE_PREFIX} \
    --bgpadding 1000 \
    -t 10
