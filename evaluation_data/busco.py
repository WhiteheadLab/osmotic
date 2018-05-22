import os
import os.path
from os.path import basename
import subprocess
from subprocess import Popen, PIPE
# custom Lisa module
import pandas as pd



def parse_busco_stats(busco_filename,sample):
        print(busco_filename)
        count=0
        important_lines=[10,11,12,13,14,15]
        busco_dict={}
        busco_dict[sample]=[]
        if os.stat(busco_filename).st_size != 0:
            with open(busco_filename) as buscofile :
                for line in buscofile:
                    count+=1
                    line_data=line.split()
                    print(count)
                    print(line_data)
                    if count in important_lines:
                        if count == 8:
                            line_data = line_data[0].split(":")
                            line_data = line_data[1].split("[")
                            line_data = line_data[0].split("%")
                            print(line_data)
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

def execute(data_frame,busco_dir,species,basedir):
	# construct an empty pandas dataframe to add on each assembly.csv to
        newdir=basedir+busco_dir+"/"
        files = os.listdir(newdir)
        for filename in files:
            if filename.startswith("short_summary"):
                busco_file = newdir + filename
		#run_busco(busco_dir,fixed_trinity_fasta,species)
                data=parse_busco_stats(busco_file,species)
		#data_frame=build_DataFrame(data_frame,data)
        return data_frame

data_frame=pd.DataFrame()
basedir = "/home/ljcohen/osmotic_killifish/euk_busco/"
listofdirs=os.listdir(basedir)
print(listofdirs)
for busco_dir in listofdirs:
    species = busco_dir.split(".")[0][5:]
    print(species)
    data_frame=execute(data_frame,busco_dir,species,basedir)

data_frame.to_csv("busco_scores_v3_euk.csv")
print("BUSCO stats written.")                     
