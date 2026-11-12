# scRNA-seq: build the GRCh38 reference with EGFP_EXTRA and TDTOMATO_EXTRA_IVS
# reporter transcripts added (cellranger 3.0.2 mkref). Run by the BU CReM
# bioinformatics core. Reporter features appear in Seurat as "EGFP-EXTRA" and
# "TDTOMATO-EXTRA-IVS" (Seurat converts "_" to "-").

cd /restricted/projectnb/crem-bioinfo/reference_data/REFERENCES/GRCh38_tdtomato_GFP_10X
cp /restricted/projectnb/crem-bioinfo/reference_data/REFERENCES//GRCh38//Homo_sapiens.GRCh38.92.gtf Homo_sapiens.GRCh38.93.filtered_GFP_Tomato.gtf
cp /restricted/projectnb/crem-bioinfo/reference_data/REFERENCES//GRCh38//Homo_sapiens.GRCh38.dna.primary_assembly.fa Homo_sapiens.GRCh38.dna.primary_assembly_GFP_Tomato.fa

echo -e 'TDTOMATO_EXTRA_IVS\tTDTOMATO_EXTRA_IVS_marker\texon\t1\t1746\t0.0\t+\t.\tgene_id "TDTOMATO_EXTRA_IVS"; transcript_id "TDTOMATO_EXTRA_IVS"; gene_name "TDTOMATO_EXTRA_IVS"; gene_biotype "protein_coding";' > TDTOMATO_EXTRA_IVS.gtf
echo -e 'EGFP_EXTRA\tEGFP_EXTRA_marker\texon\t1\t756\t0.0\t+\t.\tgene_id "EGFP_EXTRA"; transcript_id "EGFP_EXTRA"; gene_name "EGFP_EXTRA"; gene_biotype "protein_coding";' > EGFP_EXTRA.gtf
cat EGFP_EXTRA.gtf >> Homo_sapiens.GRCh38.93.filtered_GFP_Tomato.gtf
cat TDTOMATO_EXTRA_IVS.gtf >> Homo_sapiens.GRCh38.93.filtered_GFP_Tomato.gtf
cat EGFP_EXTRA.fa >> Homo_sapiens.GRCh38.dna.primary_assembly_GFP_Tomato.fa
cat TDTOMATO_EXTRA_IVS.fa >> Homo_sapiens.GRCh38.dna.primary_assembly_GFP_Tomato.fa

grep ">" Homo_sapiens.GRCh38.dna.primary_assembly_GFP_Tomato.fa
tail Homo_sapiens.GRCh38.93.filtered_GFP_Tomato.gtf

module load bcl2fastq/2.20
module load cellranger/3.0.2
cellranger mkref --genome=GRCh38_tdtomato_GFP_10X --fasta=Homo_sapiens.GRCh38.dna.primary_assembly_GFP_Tomato.fa --genes=Homo_sapiens.GRCh38.93.filtered_GFP_Tomato.gtf

