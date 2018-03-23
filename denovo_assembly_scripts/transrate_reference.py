import os
import os.path
from os.path import basename
# custom Lisa module
import clusterfunc_py3
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


def transrate(transratedir,transrate_out,trinity_fasta,sample,reference):
	transrate_command = """
transrate --assembly={} --reference={} \
--threads=4 \
--output={}
""".format(trinity_fasta,reference,transrate_out)
    	print transrate_command
    	commands = [transrate_command]
    	process_name = "transrate"
   	module_name_list = ""
    	filename = sample
    	clusterfunc_py3.sbatch_file(transratedir, process_name,module_name_list, filename, commands)

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

def execute(data_frame, listoffiles, assemblydir,transratedir,reference):
	#pairs_dictionary=get_pairs(listoffiles,basedir)
    	# construct an empty pandas dataframe to add on each assembly.csv to
    	for fasta in listoffiles:
            if fasta.endswith(".fasta"):
                sample = fasta.split(".")[0]
        	trinity_fasta = assemblydir + fasta
		transrate_out = transratedir + sample + "/"
        	transrate_assemblies = transrate_out + "/" + "assemblies.csv"
		if os.path.isfile(transrate_assemblies):
        	        data = parse_transrate_stats(transrate_assemblies)
                    	data_frame = build_DataFrame(data_frame, data)
        	else:  
                    	print "Running transrate..."
                  	transrate(transratedir,transrate_out,trinity_fasta,sample,reference)
	return data_frame

assemblydir = "/home/ljcohen/public_html/killifish/transcriptome_assemblies/"
transratedir = "/home/ljcohen/osmotic_transrate/kfish2rae5g_mrna.combined/"
reference = "/home/ljcohen/reference/kfish2rae5/kfish2rae5g.mrna.combined"
clusterfunc_py3.check_dir(transratedir)
listoffiles = os.listdir(assemblydir)
data_frame = pd.DataFrame()
data_frame = execute(data_frame,listoffiles, assemblydir, transratedir,reference)
#data_frame.to_csv("transrate_scores.csv")
