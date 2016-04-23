import fnmatch
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

#1. Get data from spreadsheet

def get_files(assemblydirs,basedir):
        files_dictionary={}
	for genus_species in assemblydirs:
		genus_species_dir=assemblydir+genus_species+"/"
		listoffiles=os.listdir(genus_species_dir)
        	for filename in listoffiles:
			#print filename
                	if filename.endswith("left.fq"):
				full_filename=genus_species_dir+filename
                        	fields=filename.split(".")
                        	sample_name=fields[0]
				if sample_name in files_dictionary.keys():
                                	files_dictionary[sample_name].append(full_filename)
                        	else:
                                	files_dictionary[sample_name]=[full_filename]
        return files_dictionary

def get_cat_command(genus_species_dir,genus_species):
	cat_command="""
cat {}*.1 > {}{}.fix.left.fq
cat {}*.2 > {}{}.fix.right.fq
gunzip -c {}orphans.keep.abundfilt.fq.gz >> {}{}.left.fq
""".format(genus_species_dir,genus_species_dir,genus_species,genus_species_dir,genus_species_dir,genus_species,genus_species_dir,genus_species_dir,genus_species)
	return cat_command

def get_sourmash_command(full_filename,genus_species):
	sourmash_command="""
head -400000000 {} > /home/ljcohen/sourmash_temp/{}.head
""".format(full_filename,genus_species)
# /home/ubuntu/sourmash/sourmash compute --protein -k 18,21 -
	print sourmash_command
	#s=subprocess.Popen(sourmash_command,shell=True)
        #s.wait()

        	
def execute(assemblydirs,assemblydir):
        #files_dictionary=get_files(assemblydirs,assemblydir)
	for genus_species in assemblydirs:
		genus_species_dir=assemblydir+genus_species+"/"
		cat_command=get_cat_command(genus_species_dir,genus_species)			
		command=[cat_command]
		module_load_list=[""]
        	process_name="cat"
        	clusterfunc.sbatch_file(genus_species_dir,process_name,module_load_list,genus_species,command)



			#s=subprocess.Popen(sourmash_command,shell=True)
    			#s.wait()
			#if os.path.isfile("*.sig"):
			#	print os.path.listdir(trimdir)
			#else:
			#	print "sourmash not run yet"
				

assemblydir="/home/ljcohen/osmotic_assemblies_completed/"
assemblydirs=os.listdir(assemblydir)
execute(assemblydirs,assemblydir)
