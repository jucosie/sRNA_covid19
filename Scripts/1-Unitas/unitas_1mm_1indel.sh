#!/bin/bash
#SBATCH --job-name=unitas_1mm_1indel								# Job name.
#SBATCH --output=unitas_1mm_1indel.log                                       		  	# Standard output and error log.
#SBATCH --qos medium								# Partition (queue)
#SBATCH --ntasks=1									# Run on one mode. 
											# Don't change unless you know what you are doing.
#SBATCH--cpus-per-task=1								# Number of tasks = cpus. 
											# It depends on the number of process of your parallelization
#SBATCH --time=1-00:00:00								# Time limit days-hrs:min:sec.
#SBATCH --mem=60gb									# Job memory request.

input_dir="path/to/your/input/dir"
results_unitas="path/to/your/output/dir"

cd ${results_unitas}

#####################################################
# This program performs the adapter trimming, size  #
# selection and annotation of raw sRNA libraries    #
#####################################################

srun perl /storage/ncRNA/Softwares/unitas_1.8.0.pl -input ${input_dir} -species homo_sapiens -trim AGATCGGAAGAG -skip_dust off -trim_minlength 12 -trim_maxlength 34 -keep_temp -species_miR_only -na_for_phasi -tail 0 -insdel 1
