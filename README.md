# BRC_Organoid_Joana
Transcriptomic changes in organoid small intestinal epithelial cells upon interaction with lymphocytes

**DB ID** Project.id=6, Data.id=11

## Scripts
1. header.R
  * header file with global variables and settings
2. 01_createDbEntry.R 
  * create database entries for the samples  
3. 02_fastQC.R  
  * QA checks on the FASTQ files using the class CFastqQuality  
4. 03_trimmomatic_array_job.R  
  * writes a batch array job for hpc to use trimmomatic for read trimming  
5. 04_fastQC_trimmed.R  
  * QA checks on the trimmed FASTQ files using the class CFastqQuality  
6. 05_hisat2_array_job.R  
  * create a parameter file and shell script to run array job on hpc  
7. 06_samtools_array_job.R  
  * create a parameter file and shell script to run array job on hpc  
8. 07_bam_files_qa.R  
  * quality checks on the bam files  
9. 08_counts_from_bams.R  
  * generate count tables for transcripts from bam files  
10. 09_clustering_counts_matrix.R  
  * cluster the samples based on count table  
11. 10_de_analysis.R  
  * DE analysis for the count data using DESeq2  
12. 11_de_analysis_MCMC.R  
  * DE analysis for the count data using a Bayesian approach by writing a log posterior function and optimizing the function using LearnBayes::laplace. The second half of the script has the option to use Stan to generate MCMC samples from the posterior.  
  * **nb_glm.stan** Stan script, similar to the log posterior function. Can be used for genes that do not converge using optimizer.  

  
