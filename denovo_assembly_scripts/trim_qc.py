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

def parse_filename(filename):
	listoffilestomerge=[]
	fields=filename.split("_")
	genus=fields[0]
	species=fields[1]
	population=fields[2]
	treatment=fields[3]
	sample=fields[4]
	read=fields[5]
	extension=fields[6]
	sample=(genus,species,population,treatment,sample)
	return sample

def get_trimmomatic(trimdir,file1,file2,sample):
	j="""#!/bin/bash
java -Xmx10g -jar /home/nreid/bin/trimmomatic-0.33/trimmomatic-0.33.jar PE \\
-baseout {}{}.trim.fq \\
{} {} \\
ILLUMINACLIP:/home/nreid/rnaseq/NEBnextAdapt.fa:2:40:15 \\
SLIDINGWINDOW:4:2 \\
LEADING:2 \\
TRAILING:2 \\
MINLEN:25 &> {}trim.{}.log
""".format(trimdir,sample,file1,file2,trimdir,sample)
	return j

def fastqc_report(fastqcdir,fastq_file,sample_name):
        fastqc_string="fastqc -o "+fastqcdir+" "+fastq_file
	process_string=[fastqc_string]
	process_name="fastqc"
	module_load_list=["fastqc/0.10.1"]
	clusterfunc.sbatch_file(basedir,process_name,module_load_list,sample_name,process_string)

def get_trimmed_list(listoffiles):
	fastq_file_list=[]
	for filename in listoffiles:
        	if filename.endswith(".TruSeq_1U.fq"):
			fastq_file_list.append(fastqcdir+filename)
        	elif filename.endswith(".TruSeq_2U.fq"):
                	fastq_file_list.append(fastqcdir+filename)

def run_trimmomatic(trimdir,file1,file2,sample):
	trimmomatic_string=get_trimmomatic(trimdir,file1,file2,sample)
	process_string=[trimmomatic_string]
	module_load_list=""
	process_name="trim"
	#clusterfunc.sbatch_file(trimdir,process_name,module_load_list,sample,process_string)

def get_pairs(listoffiles,basedir):
	pairs_dictionary={}
	for basefilename in listoffiles:
		if basefilename.endswith(".fastq.gz"):
                	filename=basedir+basefilename
                	fields=basefilename.split("_")
                	sample_name_info=fields[:-1]
			sample_name="_".join(sample_name_info)
			if sample_name in pairs_dictionary.keys():
				pairs_dictionary[sample_name].append(basedir+basefilename)
			else:
				pairs_dictionary[sample_name]=[basedir+basefilename]
	return pairs_dictionary
	
def execute(listoffiles,basedir,fastqcdir,trimdir):
	pairs_dictionary=get_pairs(listoffiles,basedir)
	for sample in pairs_dictionary.keys():
		print sample
		fileslist=sorted(pairs_dictionary[sample])
		print fileslist
		for filename in fileslist:
			fastqc_report(fastqcdir,filename,sample)
	        R1=fileslist[0]
		R2=fileslist[1]
		run_trimmomatic(trimdir,R1,R2,sample)	


basedir="/home/ljcohen/osmotic_combined/"
trimdir="/home/ljcohen/osmotic_trim/"
fastqcdir="/home/ljcohen/osmotic_trim/fastqc/"
listoffiles=os.listdir(basedir)
execute(listoffiles,basedir,fastqcdir,trimdir)
