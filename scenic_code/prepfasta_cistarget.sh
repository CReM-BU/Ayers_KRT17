
REGION_BED="/Users/ljayers/Desktop/scenicoutput/ATAC_06_24_25/region_sets/ArchR_allPeaks.bed"
GENOME_FASTA="/Users/ljayers/Desktop/scenic/reference/hg38.fa"
CHROMSIZES="/Users/ljayers/Desktop/scenic/reference/hg38.chrom.sizes"
DATABASE_PREFIX="ATAC_062425"
SCRIPT_DIR="/Users/ljayers/Desktop/scenic/create_cisTarget_databases"

${SCRIPT_DIR}/create_fasta_with_padded_bg_from_bed.sh \
        ${GENOME_FASTA} \
        ${CHROMSIZES} \
        ${REGION_BED} \
        ATAC_062425_padded.fa \
        1000 \
        yes