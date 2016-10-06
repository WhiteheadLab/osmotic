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
import pandas as pd

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

def fix_fasta(trinity_fasta, trinity_dir, sample):
        # os.chdir(trinity_dir)
    trinity_out = trinity_dir + sample + ".Trinity.fixed.fa"
    fix = """
sed 's_|_-_g' {} > {}
""".format(trinity_fasta, trinity_out)
    # s=subprocess.Popen(fix,shell=True)
    print fix
    # s.wait()
    # os.chdir("/mnt/home/ljcohen/MMETSP/")
    return trinity_out


def transrate(transratedir,transrate_out,trinity_fasta,sample,left,right):
	transrate_command = """
transrate --assembly={} --threads=4 \
--left={} \
--right={} \
--output={}
""".format(trinity_fasta,left,right,transrate_out)
    	print transrate_command
    	commands = [transrate_command]
    	process_name = "transrate"
   	module_name_list = ""
    	filename = sample
    	clusterfunc.sbatch_file(transratedir, process_name,module_name_list, filename, commands)

def parse_transrate_stats(transrate_assemblies):
    print transrate_assemblies
    if os.stat(transrate_assemblies).st_size != 0:
        data = pd.DataFrame.from_csv(transrate_assemblies, header=0, sep=',')
        return data

def build_DataFrame(data_frame, transrate_data):
    # columns=["n_bases","gc","gc_skew","mean_orf_percent"]
    frames = [data_frame, transrate_data]
    data_frame = pd.concat(frames)
    return data_frame

def execute(data_frame, listoffiles, assemblydir,transratedir):
	samples = os.listdir(assemblydir)
	sample_dictionary = {}
	#pairs_dictionary=get_pairs(listoffiles,basedir)
    	# construct an empty pandas dataframe to add on each assembly.csv to
    	for sample in samples:
        	if sample.startswith("F_heteroclitus"):
			diginormdir = "/home/ljcohen/osmotic_trim_F_heteroclitus/diginorm/"
		elif sample.startswith("F_diaphanus"):
			diginormdir = "/home/ljcohen/osmotic_trim_F_diaphanus/diginorm/"
		elif sample.startswith("F_sciadicus"):
			diginormdir = "/home/ljcohen/osmotic_trim_F_sciadicus/diginorm/"
		else:
			diginormdir = "/home/ljcohen/osmotic_trim/assemblies_msu/"+sample+"/"
		if os.path.isdir(diginormdir):
			#print diginormdir
			left = diginormdir + sample +".left.fq"
			right = diginormdir + sample + ".right.fq"
			if os.path.isfile(left) and os.path.isfile(right):
				print left
				print right
			else:
				print "Does not exist.",left,right
		else:
			print "Does not exist:",diginormdir
        	trinity_fasta = assemblydir + sample + "/" + sample + ".Trinity.fixed.fa"
		transrate_out = transratedir + sample + "/"
        	transrate_assemblies = transrate_out + "/" + "assemblies.csv"
		if os.path.isfile(trinity_fasta):
        		print trinity_fasta
        	else:
                	print "Trinity failed:", trinity_fasta
		if os.path.isfile(transrate_assemblies):
        	        data = parse_transrate_stats(transrate_assemblies)
                    	data_frame = build_DataFrame(data_frame, data)
        	else:  
                    	print "Running transrate..."
                  	transrate(transratedir,transrate_out,trinity_fasta,sample,left,right)
	return data_frame

assemblydir = "/home/ljcohen/msu_assemblies_finished/"
basedir = "/home/ljcohen/osmotic_combined/"
transratedir = "/home/ljcohen/osmotic_transrate_scores/"
clusterfunc.check_dir(transratedir)
listoffiles = os.listdir(basedir)
data_frame = pd.DataFrame()
data_frame = execute(data_frame,listoffiles, assemblydir, transratedir)
#data_frame.to_csv("transrate_scores.csv")
