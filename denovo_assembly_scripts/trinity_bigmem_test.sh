#!/bin/bash -l
#SBATCH -D /home/ljcohen/slurm_out/
#SBATCH -J trinity_F_grandis
#SBATCH -A millermrgrp
#SBATCH -p bigmeml
#SBATCH -t 96:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32
#SBATCH --mem=260000

module load rsem/1.2.23
module load trinity/2.2.0

PROJECTDIR=/home/ljcohen/assemblies/
SPECIES=F_grandis

mkdir /scratch/$SLURM_JOBID
cd /scratch/$SLURM_JOBID
cp $PROJECTDIR/$SPECIES/*.fq .

Trinity --left ./$SPECIES.left.fq --right ./$SPECIES.right.fq --output /scratch/$SLURM_JOBID/$SPECIES.trinity_out --full_cleanup --seqType fq --max_memory 256G --CPU 30
cp /scratch/$SLURM_JOBID/$SPECIES.trinity_out.Trinity.fasta /home/ljcohen/osmotic_assemblies_farm/
rm -rf /scratch/$SLURM_JOBID*
