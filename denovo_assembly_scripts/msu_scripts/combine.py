# combines files from same sample
# first 9 fields must be the same for F. heteroclitus
# first 8 fields must b the same for all other species

import os
import os.path
import subprocess
from subprocess import Popen,PIPE
import glob
import gzip
import clusterfunc

#dir1="/home/nreid/rnaseq/rawdata/141212_HS2A/"
#dir2="/home/nreid/rnaseq/rawdata/141203_HS4B/"
#dir3="/home/nreid/rnaseq/rawdata/141126_HS3A/"

def sort_filenames_by_genus_species_treatment(listoffiles,basedir):
	files={}
	for i in listoffiles:
		if i.endswith(".fastq.gz"):
			genus_species_sample_read=parse_filename(i)
			if genus_species_sample_read in files:
				files[genus_species_sample_read].append(i)
			else:
				files[genus_species_sample_read]=[i]
	for species in files.keys():
	#	filestomergelist=find_files_to_merge(files[species],species)
		print species
		files[species]=sorted(files[species])
		print files[species]
		print len(files[species])
	return files

def combine_files(files,basedir,combine_dir):
	clusterfunc.check_dir(combine_dir)
	for species in files.keys():
		fields=files[species][0].split("_")
		extension=fields[-1]
		parsed_extension1=extension.split(".")
		parsed_extension2=parsed_extension1[1:]
		new_extension=".".join(parsed_extension2)
		newfilename=get_newfilename(fields,new_extension)
		print species
		print files[species]
		#print newfilename
		newfilename=combine_dir+newfilename
		files_string=" ".join(files[species])
		combine_string="cat "+files_string+" > "+newfilename
		print combine_string
		workingdir=os.getcwd()
		os.chdir(basedir)
		#s=subprocess.Popen(combine_string,shell=True)
		#s.wait()
		os.chdir(workingdir)
		
def parse_extension(extension):
	parsed_extension1=extension.split(".")
	parsed_extension2=parsed_extension1[1:]
	filenumber=parsed_extension1[0]
	trim=parsed_extension1[1]
	new_extension=".".join(parsed_extension2)
        return filenumber,trim

def find_files_to_merge(fileslist,species):
	filestomergelist=[]
        for filename in fileslist:
		fields=filename.split("_")
	return filestomergelist

def get_newfilename(fields,new_extension):
	newfilename1="_".join(fields[0:-1])
	newfilename=newfilename1+"_"+new_extension
	return newfilename

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
		filenumber,trim=parse_extension(extension)
		genus_species_sample_read=(genus,species,population,treatment,sample,flowcell,lane,read,trim)
		return genus_species_sample_read
	else:	
		treatment=fields[3]
		sample=fields[4]
		barcode=fields[5]
		lane=fields[6]
		read=fields[7]
		extension=fields[8]
		filenumber,trim=parse_extension(extension)
		population="NA"
		genus_species_sample_read=(genus,species,population,treatment,sample,flowcell,lane,read,trim)
		return genus_species_sample_read

# how to figure out which files need to be combined?
# the last field will be different
# all the rest of the fields will be the same
# need to collect this information somehow for all files

# put all files to be merged in a list:

#basedir="/home/ljcohen/osmotic_trim/Ph30/"
basedir="/home/ljcohen/osmotic_trim/Ph40/"
listoffiles=os.listdir(basedir)
#print listoffiles
files=sort_filenames_by_genus_species_treatment(listoffiles,basedir)
combine_dir="/home/ljcohen/osmotic_trim/Ph40/combined/"
combine_files(files,basedir,combinedir)
