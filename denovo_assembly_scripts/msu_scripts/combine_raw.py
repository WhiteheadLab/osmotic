import os
import os.path
import subprocess
from subprocess import Popen,PIPE
import glob
import gzip
import clusterfunc

def parse_filename(fastq_filename):
# example filename:
# AWJRDD001_F_similis_BW_1_TCCGCGAA-ATAGAGGC_L002_R1_001.fastq.gz
# 0=flowcell#
# 1=Genus
# 2=spcies
# 3=treatment
# 4=sample#
# 5=barcode
# 6=lane
# 7=read
## this can be parsed by .split(".")
# 8=file#.fastq.gz

# F. heteroclitus:
# AWJRDD002_F_heteroclitus_MDPL_transfer_2_ATTCAGAA-ATAGAGGC_L004_R1_002.fastq.gz 
# 0=flowcell#
# 1=Genus
# 2=species
# 3=population
# 4=treatment
# 5=sample#
# 6=barcode
# 7=lane
# 8=read
## this can be parsed by .split(".")
# 9=file#.fastq.gz
	listoffilestomerge=[]
	fields=fastq_filename.split("_")
	flowcell=fields[0]
	genus=fields[1]
	species=fields[2]
	genus_species=genus+"_"+species
	if species=="heteroclitus":
		population=fields[3]
		treatment=fields[4]
		sample=fields[5]
		barcode=fields[6]
		lane=fields[7]
		read=fields[8]
		extension=fields[9]
		sample=(genus,species,population,treatment,sample,read)
		return sample
	else:	
		treatment=fields[3]
		sample=fields[4]
		barcode=fields[5]
		lane=fields[6]
		read=fields[7]
		extension=fields[8]
		population="NA"
		sample=(genus,species,population,treatment,sample,read)
		return sample

def get_newfilename(sample,out_dir):
	newfilename_base="_".join(sample)
	newfilename=out_dir+newfilename_base+".fastq.gz"
	return newfilename

def collect_files(data_dir):
	filename_dictionary={}
	for root,dirs,files in os.walk(data_dir):
		for filename in files:
			if filename.endswith(".fastq.gz"):
				sample=parse_filename(filename)
				if sample in filename_dictionary.keys():
					filename_dictionary[sample].append(root+"/"+filename)
				else:
					filename_dictionary[sample]=[root+"/"+filename]
	return filename_dictionary 	

def get_file_string(fileslist):
	sorted_fileslist=sorted(fileslist)
	files_string=" ".join(sorted_fileslist)
	return files_string

def combine_files(filename_dictionary,output_dir):
	for sample in filename_dictionary.keys():
		newfilename=get_newfilename(sample,output_dir)
		#print newfilename
		files_string=get_file_string(filename_dictionary[sample])
		combine_string="cat "+files_string+" >> "+newfilename
		print combine_string
		combine_string_submit=[combine_string]
		process_name="combine"
		module_name_list=""
		filename="_".join(sample)
		#clusterfunc.sbatch_file(output_dir,process_name,module_name_list,filename,combine_string_submit)

output_dir="/home/ljcohen/osmotic_combined/"
data_dir="/home/nreid/rnaseq/rawdata/"
filename_dictionary=collect_files(data_dir)
for sample in filename_dictionary.keys():
	print sample
	print "This is the number of files to combine:",len(filename_dictionary[sample])
	print sorted(filename_dictionary[sample])

#combine_files(filename_dictionary,output_dir)
