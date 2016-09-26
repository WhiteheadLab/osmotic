import os
import os.path
from os.path import basename
import subprocess
from subprocess import Popen, PIPE
import clusterfunc

def get_files(dirname):
	listoffiles=os.listdir(dirname)
	raw_stats="read_stats.txt"
	with open(raw_stats,"w") as statsfile:
		for filename in listoffiles:
			if filename.endswith("fastq.gz"):
				print filename
				statsfile.write(filename+"\t")
				command="zcat "+dirname+filename+" | wc -l"
				s=subprocess.Popen(command,stdout=subprocess.PIPE,shell=True)
				s.wait()
				print command
				out=s.communicate()[0]
				statsfile.write(out+"\n")
			else:
				print "File not found:",filename
		
	print "File written:",raw_stats

def get_trim_files(dirname):
        listoffiles=os.listdir(dirname)
	trim_stats="trim_stats.txt"
	with open(trim_stats,"w") as statsfile:
        	for filename in listoffiles:
			if filename.endswith("P.fq"):
                		print filename
                		statsfile.write(filename+"\t")
				command="cat "+dirname+filename+" | wc -l"
                		s=subprocess.Popen(command,stdout=subprocess.PIPE,shell=True)
                		s.wait()
                		out=s.communicate()[0]
				statsfile.write(out+"\n")
	print "File written:",trim_stats
	
def get_genus_species_reads(assembly_dir):
	diginorm_stats="diginorm_stats.txt"
	with open(diginorm_stats,"w") as statsfile:
		genus_species_dirs=os.listdir(assembly_dir)
		for genus_species in genus_species_dirs:
			genus_species_dir=assembly_dir+genus_species+"/"
			listoffiles=os.listdir(genus_species_dir)
			for filename in listoffiles:
				if filename.endswith("readcount"):
					print filename
					statsfile.write(filename+"\t")
					stats_filename=genus_species_dir+filename
					with open(stats_filename,"rU") as inputfile:
						line=next(inputfile).split(" ")
						reads=line[2].rstrip()
						statsfile.write(reads+"\n")
	print "File written:",diginorm_stats 

assembly_dir="/home/ljcohen/osmotic_trim/assemblies/"
combined_dir="/home/ljcohen/osmotic_combined/"
trim_dir="/home/ljcohen/osmotic_trim/"
get_files(combined_dir)
#get_trim_files(trim_dir)
#get_genus_species_reads(assembly_dir)
