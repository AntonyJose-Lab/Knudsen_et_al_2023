#!/bin/bash
#####
# Run this script after '1_fastq_to_sorted_bam' and '2_sorted_bam_to_mutated_genes.sh'. This script takes mutations in multiple strains and identifies potential 'complementation groups' that have different mutations in the same gene.
#####
screen_name="my_screen"
path_to_project="/path/to/project/"
cd ${path_to_project}/analyses/${screen_name}_variants/${screen_name}_annotated_vcf/new_mutations
mkdir insilico_comp

### Collects together all genes that occur more than once based on frequency.
cut -f 11 *.vcf | sort | uniq -c | awk '{ if($1 > 1) {print $2} }' > insilico_comp/multi_freq_genes.txt

### For each gene identified as being mutated in more than one strain, go through and check all lines of all strains below and print what is needed. If there are any mutants that seem to be outliers for the number of mutations they have (e.g., PCR introduction of mutations during library preparateion), they could be excluded and a new set of strains can be considered after modifying the for i in $(cat <filename>) appropriately.
for i in $(cat insilico_comp/multi_freq_genes.txt)
do
	touch insilico_comp/${i}.txt
	touch insilico_comp/temp1
	touch insilico_comp/temp2
	for j in $(cat ${path_to_project}/analyses/strain_names.txt)
	do
	number_mut=$(wc -l ${j}_new_mutations.vcf | awk '{print $1}')
		for (( k=1; k<=$number_mut; k++ ))
		do
		lhs=$(awk -F "|" -v n=$k 'NR == n {print $4}' ${j}_new_mutations.vcf)
		if [ "$lhs" == "$i" ]; then
		awk -v n=$j 'BEGIN {print n}' >> insilico_comp/temp1
		## This awk prints out gene name, DNA change, and protein change to a file.
		awk -F "|" -v n=$k 'NR == n {print $4,$10,$11}' ${j}_new_mutations.vcf >> insilico_comp/temp2
		cat insilico_comp/temp1
		cat insilico_comp/temp2
		fi
		done
	done
	num_strains=$(sort insilico_comp/temp1 | uniq | wc -l | awk '{print $1}')

	if (( $num_strains > 1 )); then
		paste insilico_comp/temp1 insilico_comp/temp2 >> insilico_comp/${i}.txt; else
		rm insilico_comp/${i}.txt
	fi
	> insilico_comp/temp1
	> insilico_comp/temp2
done
rm insilico_comp/temp*
rm insilico_comp/multi_freq_genes.txt
cd insilico_comp
cat * > ../${screen_name}_in_silico_complementation_groups.txt
cd ..
rm -r insilico_comp
