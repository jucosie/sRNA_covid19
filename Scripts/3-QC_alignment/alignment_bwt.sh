#!/bin/bash
#SBATCH --job-name=alignment_bwt	#Job name to show with squeue
#SBATCH --output=alignment_bwt_%j.out	#Output file	
#SBATCH --ntasks=20		#Maximum number of cores to use 
#SBATCH --time=1:00:00		#Time limit to execute the job (60 minutes)
#SBATCH --mem-per-cpu=1G	#Required Memory per core
#SBATCH --cpus-per-task=10       #CPUs assigned per task.
#SBATCH --qos=short             #QoS: short,medium,long,long-mem


input_dir="/path/to/libraries/04-homo_sapiens_libraries_clean" #Use the clean output of the in-house script
output_dir="/path/to/output"
human_index="/path/to/human_index/"
sars_cov2_index="/path/to/sarscov2_index/"

cd ${input_dir}
files=$(ls)

## Map against human reference
for file in ${files}; do
  name=$(echo $file | sed 's/filtered_cleaned_homo_sapiens_//g' | sed 's/.fa//g')
  >${output_dir}/human_${name}.log
  srun -n 1 -N 1 -c $SLURM_CPUS_PER_TASK -Q --exclusive bowtie ${human_index}/GRCh38_index -n 1 -l 10 -k 1 --best -f ${file} -S ${output_dir}/human_${name}.sam --un ${output_dir}/human_unal_${name}.fasta --threads $SLURM_CPUS_PER_TASK >> ${output_dir}/human_${name}.log 2>&1 &
done
wait

cd ${output_dir}
files_2=$(ls | grep "unal")

## Map against sars-cov-2 reference
for f in ${files_2}; do
  >${output_dir}/sars_${name}.log
  name=$(echo $f | sed 's/human_unal_//g' | sed 's/.fasta//g')
  srun -n 1 -N 1 -c $SLURM_CPUS_PER_TASK -Q --exclusive bowtie ${sars_cov2_index}/sars_cov2_index -n 1 -l 10 -k 1 --best -f ${f} -S ${output_dir}/sars_${name}.sam  --threads $SLURM_CPUS_PER_TASK >> ${output_dir}/sars_${name}.log 2>&1 &
done
wait

exit 0
