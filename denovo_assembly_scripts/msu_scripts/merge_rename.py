import os
import os.path
from os.path import basename
import subprocess
from subprocess import Popen, PIPE
import clusterfunc

def get_files(assembly_dir):
	abundfilt_files={}
	fileslist=os.listfiles(assembly_dir)
	for i in fileslist:
		if i.endswith(".abundfilt")
			file_info=i.split("_")
			genus_species=file_info[1]+"_"+file_info[2]
			
			if genus_species in abundfilt_files.keys():
				



assembly_dir="/home/ljcohen/osmotic_assemblies"
