import os
import os.path
from os.path import basename
import subprocess
from subprocess import Popen, PIPE
import clusterfunc

def parse_filename(filename):
        fields=filename.split("_")
        genus=fields[0]
        species=fields[1]
        population=fields[2]
	species_population=species+population
        treatment=fields[3]
        sample=fields[4]
        read_extension=fields[5]
        new_filename="_".join([genus,species_population,treatment,sample,read_extension])
	return new_filename

def get_rename(fastq_files_dir,filename,newfilename):
	rename_command="""
mv {}{} {}{}
""".format(fastq_files_dir,filename,fastq_files_dir,newfilename)
	return rename_command

def rename_F_heteroclitus(fastq_files_dir):
	fileslist=os.listdir(fastq_files_dir)
	for filename in fileslist:
		if filename.startswith("F_heteroclitus"):
			# combine genus_speciespopulation
			new_filename=parse_filename(filename)
 			print new_filename
			rename_command=get_rename(fastq_files_dir,filename,new_filename)
			print rename_command
			#s=subprocess.Popen(rename_command,shell=True)
			#s.wait()

fastq_files_dir="/home/ljcohen/osmotic_combined/"
rename_F_heteroclitus(fastq_files_dir)
