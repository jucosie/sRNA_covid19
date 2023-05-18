###################################################################
# Unitas returns one annotation file per sample, this script      #
# unifies the results in one consensus file. Also, it reformats   #
# some of the ncRNA naming to make it compatible with tabular     # 
# format and selects the first annotation result for every        # 
# sequence                                                        #
###################################################################

unitas_results="/path/to/unitas_results"
unitas_anno_simp="/path/to/annotation_results/unitas_simplified_annotation"
unitas_noann="/path/to/annotation_results/unitas_no-anno"
scripts_dir="/path/to/scripts/Annotation"

cd $unitas_results
unitas_dirs=$(ls)

##Simplify the annotation given by unitas
for dir in $unitas_dirs; do
	sample_name=$(echo $dir | sed 's/UNITAS_05-01-2023_//g' | sed 's/.fq_#1//g') #replace with the date of your results
	cat ${dir}/unitas.full_annotation_matrix.txt | sed 's/\tfrom\t/|/g' | sed 's/\tfrom/|/g'| sed 's/ from /|/g' | tail -n +2 > ${unitas_anno_simp}/unitas_${sample_name}_temp.txt
	awk '{print $1"\t"$3"\t"$4}' ${unitas_anno_simp}/unitas_${sample_name}_temp.txt > ${unitas_anno_simp}/unitas_${sample_name}.txt
	rm ${unitas_anno_simp}/unitas_${sample_name}_temp.txt
done

cd $unitas_anno_simp

cat unitas_* > unitas_concat_anno.txt
cat unitas_concat_anno.txt | sort -k 1 | uniq > unitas_complete_anno.txt

###################################################################
# For some reason we don't know unitas looks for piRNAs in        #
# the already written no-annotation file, thus the file           #
# contains piRNA sequences we have to remove for further analyses #
###################################################################

cd ${unitas_anno_simp}
cat unitas_complete_anno.txt | grep "piR-cluster" > unitas_piRNA_anno.txt
cat unitas_piRNA_anno.txt | awk '{print $1}' > unitas_piRNA_list.txt

cd ${unitas_results}

for dir in $(ls); do
	sample_name=$(echo $dir | sed 's/UNITAS_05-01-2023_//g' | sed 's/.fq_#1//g') #replace with the date of your results
	python3 ${scripts_dir}/filter_fasta_by_list_of_seqs.py ${dir}/fasta/unitas.no-annotation.fas ${unitas_anno_simp}/unitas_piRNA_list.txt > ${unitas_noann}/${sample_name}.no-annotation.fas
done
