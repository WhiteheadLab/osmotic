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


def get_data(thefile):
    count=0
    url_data={}
    with open(thefile,"rU") as inputfile:
        headerline=next(inputfile).split(',')
        #print headerline        
        position_name=headerline.index("ScientificName")
        position_reads=headerline.index("Run")
        position_ftp=headerline.index("download_path")
        for line in inputfile:
            line_data=line.split(',')
            name="_".join(line_data[position_name].split())
            read_type=line_data[position_reads]
            ftp=line_data[position_ftp]
            name_read_tuple=(name,read_type)
            if name_read_tuple in url_data.keys():
                if ftp in url_data[name_read_tuple]:
                    print "url already exists:", ftp
                else:
                    url_data[name_read_tuple].append(ftp)
            else:
                url_data[name_read_tuple] = [ftp]
        return url_data

def fix_fasta(trinity_fasta,trinity_dir,sample):
        trinity_out=trinity_dir+sample+".Trinity.fixed.fa"
        fix="""
sed 's_|_-_g' {} > {}
""".format(trinity_fasta,trinity_out)
        #s=subprocess.Popen(fix,shell=True)
        #print fix
        #s.wait()
        return trinity_out

def run_busco(busco_dir,trinity_fasta,sample):
	busco_command="""
busco -m trans -in {} \
--cpu 16 -l /mnt/research/ged/lisa/busco/metazoa -o {}.metazoa
""".format(trinity_fasta,sample)
	print busco_command
	commands = [busco_command]
        process_name = "busco"
        module_name_list = ""
        filename = sample
        clusterfunc.qsub_file(busco_dir,process_name,module_name_list,filename,commands) 	

def parse_busco_stats(busco_filename,sample):
        print busco_filename
	count=0
	important_lines=[7,8,9,10,11,12]
	busco_dict={}
	busco_dict[sample]=[]
	if os.stat(busco_filename).st_size != 0:
        	with open(busco_filename) as buscofile :
			for line in buscofile:
				count+=1
				line_data=line.split()
				print line_data
				if count in important_lines:
					busco_dict[sample].append(int(line_data[0]))
	busco_data=pd.DataFrame.from_dict(busco_dict,orient='index')
	busco_data.columns=["Complete","Complete/Single-Copy","Complete/Duplicated","Fragmented","Missing","Total"]
	busco_data['Complete_BUSCO_perc']=busco_data['Complete']/busco_data['Total']
	return busco_data

def build_DataFrame(data_frame,transrate_data):
        #columns=["sample","Complete","Fragmented","Missing","Total"]
	frames=[data_frame,transrate_data]
	data_frame=pd.concat(frames)
	return data_frame

def execute(data_frame,species,basedir):
	# construct an empty pandas dataframe to add on each assembly.csv to
	newdir=basedir+species+"/"
	trinitydir=newdir+"trinity_out/"
	busco_dir=newdir+"busco/"
	clusterfunc.check_dir(busco_dir)
	trinity_fasta=trinitydir+"Trinity.fasta"
	busco_file=busco_dir+"qsub_files/"+"run_"+species+".metazoa/short_summary_"+species+".metazoa"
	if os.path.isfile(trinity_fasta):
		fixed_trinity_fasta=fix_fasta(trinity_fasta,trinitydir,species)
		#run_busco(busco_dir,fixed_trinity_fasta,species)
		data=parse_busco_stats(busco_file,species)
		data_frame=build_DataFrame(data_frame,data)
	else:
		print "Trinity failed:",trinity_fasta
	return data_frame

basedirs = ["/mnt/research/ged/lisa/osmotic_killifish/assemblies_msu/",
"/mnt/research/ged/lisa/osmotic_killifish/fix_assemblies_msu/"]	

data_frame=pd.DataFrame()

for basedir in basedirs:
        listofdirs=os.listdir(basedir)
	print listofdirs
	for species in listofdirs:
        	data_frame=execute(data_frame,species,basedir)

data_frame.to_csv("busco_scores.csv")
print "BUSCO stats written."                     
