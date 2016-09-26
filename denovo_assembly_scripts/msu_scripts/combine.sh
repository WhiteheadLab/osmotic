#!/bin/bash -l
#SBATCH -J combine_reads
#SBATCH -p med

python combine.py
