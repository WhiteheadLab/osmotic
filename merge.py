# combines files from same sample
# first 9 fields must be the same for F. heteroclitus
# first 8 fields must b the same for all other species

import os
import os.path
import subprocess
from subprocess import Popen,PIPE
import glob
import gzip

# from Noah:

#in R, in working directory containing bam files:

##read in table containing file sizes in bytes

#read.table(pipe("ls -l *bam"))

## col 5 is bytes. figure out which file size is smallest
#which.min(fs[,5])

##get number of lines in smallest file (in this case, file 43, 7823469 reads)

#system("samtools view F_heteroclitus_MDPL_BW_3.bam | wc -l")

##get ratio of bytes per line (85.40499)
#fs[,43]/7823469

##make histogram

#hist(log(fs[,5]/85.40499,base=10),breaks=20)


#def get_basecall_stats():
#	sample_data={}
#	header=[]
#	basecall_stats_filename1="/home/ljcohen/osmotic/BC5NU8ACXX_Demultiplex_Stats.htm"
#        basecall_stats_filename2="/home/ljcohen/osmotic/C5VPYACXX_Demultiplex_Stats.htm"
#	data=open(basecall_stats_filename1).read()
#	root=PARSER.fromstring(data)
#	for i in root.getiterator():
#		if i.tag=="td":
			
			
def get_num_reads():
	dir1="/home/nreid/rnaseq/rawdata/141212_HS2A/"
	dir2="/home/nreid/rnaseq/rawdata/141203_HS4B/"
	dir3="/home/nreid/rnaseq/rawdata/141126_HS3A/"
	listofdirs=[dir1,dir2,dir3]
	reads_data={}
	for dirname in listofdirs:
		listofsamples=os.listdir(dirname)
		for samplename in listofsamples:
			if samplename.endswith(".zip"):
				print "not a sample name"
			elif samplename.endswith(".tar"):
				print "not a sample name"
			else:
				newdir=dirname+samplename+"/"
				listofsamplefiles=os.listdir(newdir)
				for filename in listofsamplefiles:
					if filename.endswith(".fastq.gz"):
						genus_species_sample_read=parse_filename(filename)
						print filename
						num_lines=sum(1 for line in gzip.open(newdir+filename))
						num_reads=int(num_lines)//4
						print num_reads
						num_reads_file=(num_reads,filename)
						if genus_species_sample_read in reads_data:
							reads_data[genus_species_sample_read].append(num_reads_file)
						else:	
							reads_data[genus_species_sample_read]=[num_reads_file]
	return reads_data

def get_stats_file():
	header1=["Genus,species","treatment","sample"]
	header2=["Genus","species","population","treatment","sample","flowcell","lane","read","num_reads","file_name"]
	statsfilename="/home/ljcohen/osmotic/read_stats.txt"
	reads_data=get_num_reads()
	print reads_data
	with open(statsfilename,"w") as statsfile:
		statsfile.write("\t".join(header2))	
		statsfile.write("\n")
		for sample in reads_data.keys():
			sample_list=list(sample)
			for filenum in reads_data[sample]:
				statsfile.write("\t".join(sample_list)+"\t")
				statsfile.write("\t".join(str(filenum))+"\n")
	print "Read stats written:",statsfilename
					
# want to put all files from same genus, spcies (tuple) in one directory
# files after that
def sort_filenames_by_genus_species_treatment(listoffiles):
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
		print files[species]	
		print len(files[species])

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
		genus_species_sample_read=(genus,species,population,treatment,sample,flowcell,lane,read)
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
		genus_species_sample_read=(genus,species,population,treatment,sample,flowcell,lane,read)
		return genus_species_sample_read

# how to figure out which files need to be combined?
# the last field will be different
# all the rest of the fields will be the same
# need to collect this information somehow for all files

# put all files to be merged in a list:


basedir="/home/ljcohen/osmotic_trim/Ph30/combined/interleave/"
#basedir="/home/ljcohen/osmotic_trim/Ph30/"
#basedir="/home/ljcohen/osmotic_trim/Ph40/"
listoffiles=os.listdir(basedir)
#print listoffiles
sort_filenames_by_genus_species_treatment(listoffiles)
get_stats_file()
