#!/bin/bash

#SBATCH --job-name=cfg								# Job name.
#SBATCH --output=centrifuge.log				# Standard output and error log.
#SBATCH --qos short								# Partition (queue)
#SBATCH --ntasks=1									# Run on one mode. 
											# Don't change unless you know what you are doing.
#SBATCH--cpus-per-task=25								# Number of tasks = cpus. 
											# It depends on the number of process of your parallelization
#SBATCH --time=1-00:00:00								# Time limit days-hrs:min:sec.
#SBATCH --mem=60gb									# Job memory request.


module load anaconda

output_dir="/path/to/exogenous_sRNAs_results"
unmapped_file="/path/to/exogenous_sRNAs_results/unannotated_sRNAs_unique_seqs.fa"

## CENTRIFUGE ANALYSIS
echo -e "\n--------------Running Centrifuge on $file--------------\n"
centrifuge  -x /storage/ncRNA/Projects/SARS_cov2/unified_pipeline/Additional_files/centrifuge_index/hpvc --min-hitlen 20 -f -U $unmapped_file --report-file ${output_dir}/unannotated_exogenous_sRNA_summary.tsv -S ${output_dir}/unannotated_exogenous_sRNA.out -p 50
