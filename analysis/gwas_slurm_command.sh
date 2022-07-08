#!/bin/bash
#SBATCH --job-name="gwas_SunRILs_2022"   #name of this job
#SBATCH -p short              #name of the partition (queue) you are submitting to
#SBATCH -N 1                  #number of nodes in this job
#SBATCH -n 8                 #number of cores/tasks in this job, you get all 20 physical cores with 2 threads per core with hyper-threading
#SBATCH -t 6:00:00           #time allocated for this job hours:mins:seconds
#SBATCH --mail-user=nalara@ncsu.edu   #enter your email address to receive emails
#SBATCH --mail-type=BEGIN,END,FAIL #will receive an email when job starts, ends or fails
#SBATCH -o "stdout.%j.%N"     # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH -e "stderr.%j.%N"     #optional, prints our standard error

#### Executable ####
date                          #optional, prints out timestamp at the start of the job in stdout file

module load R/4.0.3 #load R


##perl /lustre/project/genolabswheatphg/katie_scripts/cmd_alignPE_hisat.pl
##perl /lustre/project/hwwgru_mhk/pl_scripts/cmd_align_hisat.pl
Rscript --no-restore --quiet --no-save /home/nicolas.lara/project/GWAS_server_analysis.R



#cd /project/guedira_seq_map/nico/fastq_files/    #navigate to data folder
#fastq_list=("AGS2000_S1_L001_interleaved.fastq.gz" "ARGA051160-14LE31_S15_L001_interleaved.fastq.gz" "GA001138-8E36_S35_L001_interleaved.fastq.gz" "GA00190-7A14_S3_L001_interleaved.fastq.gz" "GA05450-EL52_S4_L001_interleaved.fastq.gz" "GA06493-13LE6_S20_L001_interleaved.fastq.gz" "HILLIARD_S27_L001_interleaved.fastq.gz" "LA09264C-P2_S24_L001_interleaved.fastq.gz" "NC08-23383_S6_L001_interleaved.fastq.gz" "NC8248-1_S27_L001_interleaved.fastq.gz" "SS8641_S19_L001_interleaved.fastq.gz")
#reference_genome=~"/reference/data/IWGSC/RefSeq_Assemblies/v2.1/iwgsc_refseqv2.1_assembly.fa"
#hisat2 -p 2 -x $reference_genome -U $fastq_list -S ./parent_test_SAM.sam

date                          #optional, prints out timestamp when the job ends
#End of file
