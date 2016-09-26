#!/bin/bash -l
#SBATCH -J read_stats
#SBATCH -p med

python merge.py
