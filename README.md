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


  
