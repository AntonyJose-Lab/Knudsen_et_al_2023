#!/bin/bash
#####
# Run this script after '1_fastq_to_sorted_bam'. This script takes multiple sorted bam files of sequences from multiple strains mapped to the C. elegans genome and identifies the mutations that are likely to affect protein function in each strain.
# Before running this script, (1) ensure that snpEff is stored in an appropriate directory (/path/to/snpEff) and (2) create a directory called 'genomes' under the 'data' folder and move the C. elegans fasta file (WBcel235.fa) and the samtools index file (WBcel235.fa.fai) into it.
#####
screen_name="my_screen"
# bkgd_detection_multiplier=3 This optional setting is hardcoded for >2 mutants with a mutation as being a background mutation.
path_to_project="/path/to/project/"
path_to_snpEff="/data/tools"
cd ${path_to_snpEff}/snpEff

## Loop through all strains to identify all variants and annotate them.
for i in $(cat ${path_to_project}/analyses/strain_names.txt)
	do
		mkdir ${path_to_project}/analyses/${screen_name}_variants/${i}
		mkdir ${path_to_project}/analyses/${screen_name}_variants/${i}/bcftools
		mkdir ${path_to_project}/analyses/${screen_name}_variants/${i}/snpEff
		bcftools mpileup -f ${path_to_project}/data/genomes/WBcel235.fa  -A ${path_to_project}/analyses/${screen_name}_variants/bowtie2_ce11_bam_files/${i}_sorted.bam | bcftools call -mv -Ob -o ${path_to_project}/analyses/${screen_name}_variants/${i}/bcftools/${i}_variants.bcf
		bcftools view ${path_to_project}/analyses/${screen_name}_variants/${i}/bcftools/${i}_variants.bcf -o ${path_to_project}/analyses/${screen_name}_variants/${i}/bcftools/${i}_all_variants.vcf
		java -Xmx4g -jar snpEff.jar -v -stats ${path_to_project}/analyses/${screen_name}_variants/${i}/snpEff/${i}_var.html WBcel235.99 ${path_to_project}/analyses/${screen_name}_variants/${i}/bcftools/${i}_all_variants.vcf > ${path_to_project}/analyses/${screen_name}_variants/${i}/snpEff/${i}_all_variants.ann.vcf
	done

## Create directories for storing the analyzed mutations.
cd ${path_to_project}/analyses/${screen_name}_variants/
mkdir ${screen_name}_annotated_vcf
cd ${path_to_project}/analyses/${screen_name}_variants/${screen_name}_annotated_vcf
mkdir intermediate_files
mkdir new_mutations

number_of_strains=$(wc -l ${path_to_project}/analyses/strain_names.txt | awk '{print $1}')

## Collect together mutations annotated as HIGH or MODERATE impact.
for i in $(cat ${path_to_project}/analyses/strain_names.txt)
do
	grep -e 'HIGH' -e 'MODERATE' ${path_to_project}/analyses/${screen_name}_variants/${i}/snpEff/${i}_all_variants.ann.vcf | grep -e 'GT:PL\t1/1' - > ${path_to_project}/analyses/${screen_name}_variants/${screen_name}_annotated_vcf/intermediate_files/${i}_high_moderate_variants.ann.vcf
	cut -f 1,2,4,5 ${path_to_project}/analyses/${screen_name}_variants/${screen_name}_annotated_vcf/intermediate_files/${i}_high_moderate_variants.ann.vcf > ${path_to_project}/analyses/${screen_name}_variants/${screen_name}_annotated_vcf/intermediate_files/${i}_HIGH_MODERATE.txt
	glue+="${path_to_project}/analyses/${screen_name}_variants/${screen_name}_annotated_vcf/intermediate_files/${i}_HIGH_MODERATE.txt "
done
## Create a unique list of all the mutations annotated as HIGH or MODERATE impact in all strains and count their number.
cat $glue | sort -k 2 | uniq > ${path_to_project}/analyses/${screen_name}_variants/${screen_name}_annotated_vcf/intermediate_files/AMJxxxx_all_HIGH_MODERATE.txt
number_of_unique_mutations=$(wc -l ${path_to_project}/analyses/${screen_name}_variants/${screen_name}_annotated_vcf/intermediate_files/AMJxxxx_all_HIGH_MODERATE.txt | awk '{print $1}')

touch ${path_to_project}/analyses/${screen_name}_variants/${screen_name}_annotated_vcf/intermediate_files/frequency # new file to store all mutation frequencies.
touch ${path_to_project}/analyses/${screen_name}_variants/${screen_name}_annotated_vcf/intermediate_files/mutations # new file to store all new mutations [Chr Pos Ref Mut].

## Compare mutations in each strain to complete list and keep only mutations that are unique in each strain.
for (( j=1; j<=$number_of_unique_mutations; j++ ))
do
	m=0
	lhs1=$(awk -v n=$j 'NR == n {print $1}' ${path_to_project}/analyses/${screen_name}_variants/${screen_name}_annotated_vcf/intermediate_files/AMJxxxx_all_HIGH_MODERATE.txt)
	lhs2=$(awk -v n=$j 'NR == n {print $2}' ${path_to_project}/analyses/${screen_name}_variants/${screen_name}_annotated_vcf/intermediate_files/AMJxxxx_all_HIGH_MODERATE.txt)
	for i in $(cat ${path_to_project}/analyses/strain_names.txt)
	do
		x=$(wc -l ${path_to_project}/analyses/${screen_name}_variants/${screen_name}_annotated_vcf/intermediate_files/${i}_HIGH_MODERATE.txt | awk '{print $1}')
		for (( k=1; k<=$x; k++ ))
		do
			rhs1=$(awk -v n=$k 'NR == n {print $1}' ${path_to_project}/analyses/${screen_name}_variants/${screen_name}_annotated_vcf/intermediate_files/${i}_HIGH_MODERATE.txt)
			rhs2=$(awk -v n=$k 'NR == n {print $2}' ${path_to_project}/analyses/${screen_name}_variants/${screen_name}_annotated_vcf/intermediate_files/${i}_HIGH_MODERATE.txt)
			if (( "$lhs1" == "$rhs1" && "$lhs2" == "$rhs2" )); then
			(( m++ ))
			fi
		done
	done
	awk -v n=$m 'BEGIN {print n}' >> ${path_to_project}/analyses/${screen_name}_variants/${screen_name}_annotated_vcf/intermediate_files/frequency
	# mut_x=$(($bkgd_detection_multiplier * $m))
	# if (( $mut_x < $number_of_strains )); then
	if (( $m < 2 )); then
	# write out the new mutations
	awk -v n=$j 'NR == n {print $0}' ${path_to_project}/analyses/${screen_name}_variants/${screen_name}_annotated_vcf/intermediate_files/AMJxxxx_all_HIGH_MODERATE.txt >> ${path_to_project}/analyses/${screen_name}_variants/${screen_name}_annotated_vcf/intermediate_files/mutations
	fi
done

## Add frequency to all the HIGH and MODERATE mutations.
paste ${path_to_project}/analyses/${screen_name}_variants/${screen_name}_annotated_vcf/intermediate_files/AMJxxxx_all_HIGH_MODERATE.txt ${path_to_project}/analyses/${screen_name}_variants/${screen_name}_annotated_vcf/intermediate_files/frequency > ${path_to_project}/analyses/${screen_name}_variants/${screen_name}_annotated_vcf/intermediate_files/AMJxxxx_all_HIGH_MODERATE_frequency.txt

## Print out the new mutations for each strain after subtracting background mutations that were likely present in the starting strain.
y=$(wc -l ${path_to_project}/analyses/${screen_name}_variants/${screen_name}_annotated_vcf/intermediate_files/mutations | awk '{print $1}')
	for i in $(cat ${path_to_project}/analyses/strain_names.txt)
	do
		z=$(wc -l ${path_to_project}/analyses/${screen_name}_variants/${screen_name}_annotated_vcf/intermediate_files/${i}_high_moderate_variants.ann.vcf | awk '{print $1}')
		touch ${path_to_project}/analyses/${screen_name}_variants/${screen_name}_annotated_vcf/intermediate_files/${i}_new_mutations.vcf
		for (( p=1; p<=$z; p++ ))
		do
			rhsA=$(awk -v n=${p} 'NR == n {print $1}' ${path_to_project}/analyses/${screen_name}_variants/${screen_name}_annotated_vcf/intermediate_files/${i}_high_moderate_variants.ann.vcf)
			rhsB=$(awk -v n=${p} 'NR == n {print $2}' ${path_to_project}/analyses/${screen_name}_variants/${screen_name}_annotated_vcf/intermediate_files/${i}_high_moderate_variants.ann.vcf)
			for (( l=1; l<=$y; l++ ))
			do
			lhsA=$(awk -v n=${l} 'NR == n {print $1}' ${path_to_project}/analyses/${screen_name}_variants/${screen_name}_annotated_vcf/intermediate_files/mutations)
			lhsB=$(awk -v n=${l} 'NR == n {print $2}' ${path_to_project}/analyses/${screen_name}_variants/${screen_name}_annotated_vcf/intermediate_files/mutations)
			if (( "$lhsA" == "$rhsA" && "$lhsB" == "$rhsB" )); then
			awk -v n=${p} 'NR == n' ${path_to_project}/analyses/${screen_name}_variants/${screen_name}_annotated_vcf/intermediate_files/${i}_high_moderate_variants.ann.vcf >> ${path_to_project}/analyses/${screen_name}_variants/${screen_name}_annotated_vcf/intermediate_files/${i}_new_mutations.vcf
			fi
			done
		done
	awk -F "|" '{print $4}' ${path_to_project}/analyses/${screen_name}_variants/${screen_name}_annotated_vcf/intermediate_files/${i}_new_mutations.vcf | paste ${path_to_project}/analyses/${screen_name}_variants/${screen_name}_annotated_vcf/intermediate_files/${i}_new_mutations.vcf - > ${path_to_project}/analyses/${screen_name}_variants/${screen_name}_annotated_vcf/new_mutations/${i}_new_mutations.vcf
	done

## Get list of the genes, DNA change and predicted protein change.
	cd ${path_to_project}/analyses/${project_name}/${screen_name}_variants/${screen_name}_annotated_vcf/new_mutations

	for i in $(cat ${path_to_project}/analyses/${project_name}/strain_names.txt)
	do
		number_of_mutations=$(wc -l ${path_to_project}/analyses/${project_name}/${screen_name}_variants/${screen_name}_annotated_vcf/intermediate_files/${i}_HIGH_MODERATE.txt | awk '{print $1}')
		touch ${i}_new_mutations_list.txt
		for (( j=1; j<=$number_of_mutations; j++ ))
		do
			awk -v n=${j} -F "|" 'NR == n {print $4, $10, $11}' ${i}_new_mutations.vcf >> ${i}_new_mutations_list.txt
		done
	done
