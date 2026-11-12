#!/bin/zsh
# scATAC-seq: demultiplex BCLs to FASTQ (cellranger-atac 2.1.0 mkfastq).
# Run by the BU CReM bioinformatics core on the SCC (SGE qsub).
project_name="25_06_24_lauren"
BASE="/restricted/projectnb/crem-bioinfo"
PW=$BASE/project_workspace/${project_name}/ATAC
PD=$BASE/project_data/25_06_24_lauren/ATAC/2025_06_19_KottonD_AyersL_10X_ATAC
PC=$BASE/project_code/${project_name}
task_name="mkfastq"

# loop variables
script_name="${task_name}"

# housekeeping variables
group="crem-seq"
scripts_dir=$PW/qsub_scripts
script=$scripts_dir/${script_name}.qsub
output_dir=$PW/calculations
log_dir=$PW/logs/$task_name/
log=$log_dir/${script_name}.log
err=$log_dir/${script_name}.err
mkdir -p $scripts_dir $log_dir $output_dir
echo > $log > $err # reset, because otherwise it appends

num_cores="16"
localmem="252" # current possible choices for SCC are 94, 125, 252 and 504       

# ~55 minutes
# Create the script
echo '#!/bin/zsh
  module load bcl2fastq
  date
  echo '$project_name'
  cd '$output_dir'
  time /restricted/projectnb/crem-bioinfo/tools/cellranger-atac-2.1.0/cellranger-atac mkfastq \
      --run='$PD' \
      --csv='$PD'/../2025_06_19_KottonD_AyersL_10X_ATAC_sample-sheet.csv \
      --localcores='$num_cores' \
      --localmem='$localmem'
  date
' > $script
chmod +x $script

# Submits
qsub -P $group -N $script_name -o $log -e $err -V \
  -pe omp $num_cores \
  -l mem_total="${localmem}G" \
  $script
