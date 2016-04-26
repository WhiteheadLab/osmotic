#!/bin/bash -l
#SBATCH -D /home/ljcohen/osmotic
#SBATCH -J annotate
#SBATCH -p hi
#SBATCH --mail-user=ljcohen@ucdavis.edu
#SBATCH --cpus-per-task=2
#SBATCH --mem=4GB

python annotation.py
