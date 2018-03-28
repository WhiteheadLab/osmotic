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

def fix_fasta(trinity_fasta, fixed_trinity_fasta):
	fix = """
sed 's_|_-_g' {} > {}
""".format(trinity_fasta, fixed_trinity_fasta)
    	s=subprocess.Popen(fix,shell=True)
    	print fix
    	s.wait()

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
    	#clusterfunc.qsub_file(transratedir, process_name,module_name_list, filename, commands)

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

def execute(data_frame, listoffiles, transratedir):
	samples = ["F_diaphanus","F_grandis","F_heteroclitus.MDPL","F_heteroclitus.MDPP","F_sciadicus","A_xenica","F_catanatus","F_chrysotus","F_notatus","F_notti","F_olivaceous","F_parvapinis","F_rathbuni","F_similis","F_zebrinus","L_goodei","L_parva"]
	#pairs_dictionary=get_pairs(listoffiles,basedir)
    	# construct an empty pandas dataframe to add on each assembly.csv to
    	for sample in samples:
        	if sample.startswith("F_heteroclitus"):
			diginormdir = "/mnt/research/ged/lisa/osmotic_killifish/fix_assemblies_msu/"+sample+"/"
			trinity_fasta = "/mnt/research/ged/lisa/osmotic_killifish/fix_assemblies_msu/"+sample+"/trinity_out/"+"Trinity.fasta"
			fixed_trinity_fasta = "/mnt/research/ged/lisa/osmotic_killifish/fix_assemblies_msu/"+sample+"/"+sample + ".Trinity.fixed.fa"
		elif sample.startswith("F_diaphanus"):
			diginormdir = "/mnt/research/ged/lisa/osmotic_killifish/fix_assemblies_msu/"+sample+"/"
			trinity_fasta = "/mnt/research/ged/lisa/osmotic_killifish/fix_assemblies_msu/"+sample+"/trinity_out/"+"Trinity.fasta"
			fixed_trinity_fasta = "/mnt/research/ged/lisa/osmotic_killifish/fix_assemblies_msu/"+sample+"/"+sample + ".Trinity.fixed.fa"
		elif sample.startswith("F_sciadicus"):
			diginormdir = "/mnt/research/ged/lisa/osmotic_killifish/fix_assemblies_msu/"+sample+"/"
			trinity_fasta = "/mnt/research/ged/lisa/osmotic_killifish/fix_assemblies_msu/"+sample+"/trinity_out/"+"Trinity.fasta"
			fixed_trinity_fasta = "/mnt/research/ged/lisa/osmotic_killifish/fix_assemblies_msu/"+sample+"/"+sample + ".Trinity.fixed.fa"
		else:
			diginormdir = "/mnt/research/ged/lisa/osmotic_killifish/assemblies_msu/"+sample+"/"
			trinity_fasta = "/mnt/research/ged/lisa/osmotic_killifish/assemblies_msu/"+sample+"/trinity_out/"+"Trinity.fasta"
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
			if os.path.isfile(fixed_trinity_fasta):
				transrate(transratedir,transrate_out,fixed_trinity_fasta,sample,left,right)
			else:
				fix_fasta(trinity_fasta,fixed_trinity_fasta)
	return data_frame

basedir = "/mnt/research/ged/lisa/osmotic_killifish/assemblies_msu/"
transratedir = "/mnt/research/ged/lisa/osmotic_killifish/transrate/"
clusterfunc.check_dir(transratedir)
listoffiles = os.listdir(basedir)
data_frame = pd.DataFrame()
data_frame = execute(data_frame, listoffiles, transratedir)
#data_frame.to_csv("transrate_scores.csv")
