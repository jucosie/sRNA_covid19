#!/bin/bash

#SBATCH --job-name=create_tables							# Job name.
#SBATCH --output=01-Prepare_fasta_libraries.log					# Standard output and error log.
#SBATCH --qos medium								# Partition (queue)
#SBATCH --ntasks=1									# Run on one mode. 
											# Don't change unless you know what you are doing.
#SBATCH--cpus-per-task=20								# Number of tasks = cpus. 
											# It depends on the number of process of your parallelization
#SBATCH --time=1-00:00:00								# Time limit days-hrs:min:sec.
#SBATCH --mem=150gb									# Job memory request.

# As detailed in the documentation of the Python program, this in-house script
# performs multiple tasks:
#    1. Checks for unique sequences removing sRNAs with indeterminacies (N) 
#       and keeping sRNAs within the length range of 12 to 34 (this lengths 
#       can be modified as needed).
#    2. Builds a database of the unique sequences present in all samples.
#    3. Computes counts matrix (absolute and RPM) for every sample and the union
#       of every sample.
 
# NOTE: This version of the program is made for experiments that have different times and
#       conditions so the summary_samples.tsv has to have the structure of the example in 
#       the repository. 
# NOTE2: The files fave to be named in the format "no_adapters_[time]_[condition]_[nsample].fq".
# NOTE3: The output of Unitas is a fasta file. This program needs the input to be in fastq format.
#        One can convert the fasta file to a dummie fastq file with the programs seqtk:
#        seqtk seq -F "#" ${file} > ${file}.fq

path_libraries="/path/to/libraries"
path_results="/path/to/results"
path_summary="/path/to/summary_samples"

module load python/3.11

python3 ./01-Prepare_fasta_libraries_and_count_tables_multiprocessing_v2.py \
--path-libraries ${path_libraries}/ \
--path-results ${path_results}/ \
--folder-initial-libraries clean_data \
--specie homo sapiens \
--prefix-specie HSA \
--max-digits 8 \
--threads 42 \
--path-summary ${path_summary}/summary_samples.tsv
