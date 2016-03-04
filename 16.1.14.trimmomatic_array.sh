#!/bin/bash -l
#SBATCH -J array_job
#SBATCH -o array_job_out_%A_%a.txt
#SBATCH -e array_job_err_%A_%a.txt
#SBATCH --array=1-728
#SBATCH -p med
#SBATCH --cpus-per-task=4

trimmo=java\ -jar\ /home/nreid/bin/trimmomatic-0.33/trimmomatic-0.33.jar

#fq root
fq1=$(find /home/nreid/rnaseq/ -type f -name "*.gz" | \
grep R1 | \
grep -v "ut.fastq" | \
grep -v "pt.fastq" | \
sed -n $(echo $SLURM_ARRAY_TASK_ID)p)

echo $fq1
fq2=$(echo $fq1 | sed 's/_R1_/_R2_/')
echo $fq2

ofq1=$(echo $fq1 | sed 's/fastq/pt\.fastq/' | sed 's/.*\//\/home\/ljcohen\/osmotic_trim\/Ph40\//')
ofq1u=$(echo $fq1 | sed 's/fastq/ut\.fastq/'| sed 's/.*\//\/home\/ljcohen\/osmotic_trim\/Ph40\//') 
ofq2=$(echo $fq2 | sed 's/fastq/pt\.fastq/' | sed 's/.*\//\/home\/ljcohen\/osmotic_trim\/Ph40\//')
ofq2u=$(echo $fq2 | sed 's/fastq/ut\.fastq/'| sed 's/.*\//\/home\/ljcohen\/osmotic_trim\/Ph40\//')

echo $ofq1
echo $ofq1u
echo $ofq2
echo $ofq2u

$trimmo PE -threads 4 -phred33 $fq1 $fq2 $ofq1 $ofq1u $ofq2 $ofq2u ILLUMINACLIP:NEBnextAdapt.fa:2:40:15 LEADING:2 TRAILING:2 SLIDINGWINDOW:4:2 MINLEN:25
