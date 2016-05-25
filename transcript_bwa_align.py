import os
import os.path
from os.path import basename
from urllib import urlopen
from urlparse import urlparse
import subprocess
from subprocess import Popen, PIPE
import urllib
import shutil
import glob
# custom Lisa module
import clusterfunc
# Python plotting libraries
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats, integrate


def get_assemblies(assemblydir):
	#genus_species_dirs=os.listdir(assemblydir)
	genus_species_dirs=["F_sciadicus/F_sciadicus.trinity.2/","F_diaphanus/F_diaphanus.trinity.2/"]
	for genus_species_dir in genus_species_dirs:
		genus_species=genus_species_dir.split("/")[0]
		print genus_species
		genus_species_dir=assemblydir+genus_species_dir
		assemblyfile=genus_species_dir+"Trinity.fasta"
		bam_out=genus_species_dir+genus_species+".bam"
		flagstat_out=genus_species_dir+genus_species+".flagstat.txt"
		module_load_list=["bwa/0.7.9a","samtools/1.2"]
		bwa_command=[bwa_mem(assemblyfile,bam_out,flagstat_out)]
		process_name="bwa"
		clusterfunc.sbatch_file(genus_species_dir,process_name,module_load_list,genus_species,bwa_command)

		
def bwa_index():
	bwa_index_command="""
bwa index -p killifish20130322asm /home/ljcohen/reference/killifish20130322asm.fa
""".format()
	return bwa_index_command 

def bwa_mem(assemblyfile,bam_out,flagstat_out):
	bwa_command="""
bwa mem /home/ljcohen/reference/killifish20130322asm {} | samtools view -b - > {} 
samtools flagstat {} > {}
""".format(assemblyfile,bam_out,bam_out,flagstat_out)
	return bwa_command

assemblydir="/home/ljcohen/msu_assemblies_finished/"
get_assemblies(assemblydir)
