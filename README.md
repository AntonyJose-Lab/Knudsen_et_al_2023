# Target-specific requirements for RNA interference can be explained by a single regulatory network

Daphne R. Knudsen<sup>1</sup>, Pravrutha Raman<sup>1</sup>, Farida Ettefa<sup>1</sup>, Laura De Ravin<sup>1</sup>, and Antony M. Jose<sup>1,*</sup>

<sup>1</sup>Department of Cell Biology and Molecular Genetics, University of Maryland, College Park, USA.

<sup>*</sup>Corresponding author email:  amjose@umd.edu

This paper describes some of the results from a genetic screen, where the mutated genes were identified using whole genome sequencing. The resultant data were analyzed using the uploaded shell scripts. For each strain sequenced, the fastq sequences were trimmed using cutadapt (v. 3.5), mapped to the C. elegans genome (WBcel235/ce11) using bowtie2 (v. 2.4.2), sorted using samtools (v. 1.11), and the resulting .bam file was analyzed to call variants using snpEff (v. 5.0e). The variants classified as ‘HIGH’ or ‘MODERATE’ in the .ann.vcf file for each strain that are not shared by any two or more strains were culled as new mutations caused by mutagenesis in each strain. These new mutations in each strain were compared with those of all other strains (‘in silico complementation’) using a custom script to identify sets of strains with different mutations in the same genes. 

Details for each step are commented on within the scripts. Modify paths as required before use. 

    1_fastq_to_sorted_bam.sh

    2_sorted_bam_to_mutated_genes.sh

    3_in_silico_complementation.sh

The paper also explores the response to ingested double-stranded RNA using three models of increasing complexity. One, a single-network model of protein factors with branching pathways for RNA amplification and subsequent gene silencing. Two, an equilibrium model for the dependence of mRNA and pre-mRNA on small RNAs and other RNA intermediates. Three, a dynamic model using ordinary differential equations for the dependence of mRNA and pre-mRNA on small RNAs and other RNA intermediates. Simulations of the single network and exploration of the equilibrium model were conducted in R (v. 3.6.3). Simulations of the dynamic model were conducted in Python (v. 3.8.5) and in R (v. 4.1.0).

Details for each step of the program are commented on within each program. Modify paths as required before use.

Programs in R:

    (1) For single-network model

      2022_2_9_RNAi_network_thresholds_simpler.R

    (2) For equilibrium model

      2022_6_13_RNAi_in_Celegans_linear_modified_parameter_subspace_premRNA_increase_good_kd.R

      2022_6_13_RNAi_in_Celegans_linear_modified.R

    (3) For dynamic model

      2022_7_14_RNAiDynamics_ODEs_Parameter_Analysis.R

      2022_7_27_RNAiDynamics_ODEs_Parameter_Analysis_non_SS.R

      2022_6_29_Celegans_RNAi_ODEs_RK4_method_d6.py_by_AMJ_adapted_by_DRK.R

Programs in Python:

    For dynamic model

      2022_6_29_Celegans_RNAi_ODEs_RK4_method_d6_param_sweep_non_steady_state.py

      2022_6_29_Celegans_RNAi_ODEs_RK4_method_d6_param_sweep.py

      2022_6_29_Celegans_RNAi_ODEs_RK4_method_d6.py

      2022_7_29_Celegans_RNAi_ODEs_RK4_method_d6_no_secondary_small_RNA.py

