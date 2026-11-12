#!/bin/zsh
# scATAC-seq: align and count each library (cellranger-atac 2.1.0 count)
# against refdata-cellranger-arc-GRCh38-2020-A-2.0.0 (10x stock reference;
# no GFP/tdTomato added). Produces outs/fragments.tsv.gz used by
# 03_scATAC_analysis/01_build_archr_project.R.

# Test changes
project_name="25_06_24_lauren"
BASE="/restricted/projectnb/crem-bioinfo"
PW=$BASE/project_workspace/${project_name}/ATAC
script_name="cellranger_count"

num_cores=8
localmem=128 # 16 cores * 6G = 96G is the minimum RAM recommendation (6G per core). Min

GENOME_DIR=/restricted/projectnb/crem-bioinfo/reference_data/REFERENCES/GRCh38_ATAC/refdata-cellranger-arc-GRCh38-2020-A-2.0.0
output_dir=$PW/calculations/$script_name
input_dir=/restricted/projectnb/crem-bioinfo/project_workspace/25_06_24_lauren/ATAC/calculations/AAH3NNCM5/outs/fastq_path/AAH3NNCM5

## Run cellranger

for sample_dir in ${input_dir}/* ; do
  # Set up housekeeping variables
  output_prefix=`basename ${sample_dir}`
  name=${output_prefix}
  group="crem-seq"
  scripts_dir=$PW/qsub_scripts/$script_name
  script=$scripts_dir/${name}.qsub
  log_dir=$PW/logs/$script_name
  log=$log_dir/${name}.log
  err=$log_dir/${name}.err
  mkdir -p $scripts_dir $log_dir $output_dir
  echo > $log > $err # reset, because otherwise it appends


  # Create the script
  echo '#!/bin/zsh
module load bcl2fastq/2.20

date
cd '$output_dir'
time /restricted/projectnb/crem-bioinfo/tools/cellranger-atac-2.1.0/cellranger-atac count \
  --id='$output_prefix' \
  --reference='$GENOME_DIR' \
  --fastqs='${sample_dir}' \
  --localcores='$num_cores' \
  --localmem='$localmem'
date
  ' > $script
  chmod +x $script

  # Submit
  qsub -P $group -N $name -o $log -e $err -V \
    -pe omp $num_cores \
    -l mem_total="${localmem}G" \
    $script
done
