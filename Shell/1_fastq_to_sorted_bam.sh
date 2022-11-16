#!/bin/bash
#####
# This script takes gzipped fastq files from paired-end Illumina sequencing of strains from a screen, maps them to the C. elegans genome using Bowtie2, and generates a sorted bam file.
# Before running this script, decide on a name for the screen and on the location of the project (/path/to/project/).
# Then create a directory called 'data' and store within it all gzipped fastq files named with the convention '<strain_name>_R1.fastq.gz' and '<strain_name>_R2.fastq.gz' for the paired reads. E.g., AMJ1000_R1.fastq.gz and AMJ1000_R2.fastq.gz.
# Create another folder under data called 'WBcel235' and store the bowtie2 index files wih the base name ce11 within it (i.e., ce11.rev.2.bt2, ce11.rev.1.bt2, ce11.4.bt2, ce11.3.bt2, ce11.2.bt2, ce11.1.bt2).
#####

## Change screen name as needed.
screen_name="my_screen"
## Change path_to_projects as needed.
path_to_project="/path/to/project/"
cd ${path_to_project}
## The 'path_to_project' variable stores the path to the 'data' directory
mkdir analyses
mkdir analyses/cutadapted_sequences
mkdir analyses/${screen_name}_variants/
cd analyses/${screen_name}_variants/
mkdir cutadapt_reports
mkdir bowtie2_reports
mkdir bowtie2_ce11_bam_files
## strain_names.txt below is a list of the strain names being sequenced provided as one name per line.
for i in $(cat ${path_to_project}/analyses/strain_names.txt)
do
mkdir ${path_to_project}/analyses/cutadapted_sequences/${i}_cut/
## change adapter sequences as needed.
cutadapt -j 0 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -m 20 -q 20 -o ${path_to_project}/analyses/cutadapted_sequences/${i}_cut/${i}_cut_R1.fastq.gz -p ${path_to_project}/analyses/cutadapted_sequences/${i}_cut/${i}_cut_R2.fastq.gz ${path_to_project}/data/${i}_R1.fastq.gz ${path_to_project}/data/${i}_R2.fastq.gz 1> ${path_to_project}/analyses/${screen_name}_variants/cutadapt_reports/${i}.txt
## ce11 is the base name for the Bowtie2 index files stored under the directory WBcel235.
bowtie2 -x ${path_to_project}/data/WBcel235/ce11 -1 ${path_to_project}/analyses/cutadapted_sequences/${i}_cut/${i}_cut_R1.fastq.gz -2 ${path_to_project}/analyses/cutadapted_sequences/${i}_cut/${i}_cut_R2.fastq.gz 2> ${path_to_project}/analyses/${screen_name}_variants/bowtie2_reports/${i}.txt | samtools view -bS - | samtools sort - > ${path_to_project}/analyses/${screen_name}_variants/bowtie2_ce11_bam_files/${i}_sorted.bam
done
