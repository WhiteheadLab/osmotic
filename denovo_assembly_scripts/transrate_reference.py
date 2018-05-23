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


def transrate_forward(transratedir,transrate_out,trinity_fasta,sample,reference):
	transrate_command = """
/home/ljcohen/bin/transrate-1.0.3-linux-x86_64/transrate --assembly={} --reference={} \
--threads=4 \
--output={}
""".format(trinity_fasta,reference,transrate_out)
        print(transrate_command)
    	commands = [transrate_command]
    	process_name = "transrate_forward"
   	module_name_list = ""
    	filename = sample + "_forward"
    	clusterfunc_py3.sbatch_file(transratedir, process_name,module_name_list, filename, commands)

def transrate_reverse(transratedir,transrate_out,trinity_fasta,sample,reference):
    transrate_command = """
/home/ljcohen/bin/transrate-1.0.3-linux-x86_64/transrate --assembly={} --reference={} \
--threads=4 \
--output={}
""".format(reference,trinity_fasta,transrate_out)
    print(transrate_command)
    commands = [transrate_command]
    process_name = "transrate_reverse"
    module_name_list = ""
    filename = sample + "_reverse"
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

def execute(data_frame1, data_frame2,listoffiles, assemblydir,transratedir,reference):
    #pairs_dictionary=get_pairs(listoffiles,basedir)
    # construct an empty pandas dataframe to add on each assembly.csv to
    for fasta in listoffiles:
        if fasta.endswith(".fasta"):
            print(fasta)
            sample = fasta.split(".")[0]
            trinity_fasta = assemblydir + fasta  
	    transrate_out_forward = transratedir + sample + "_trinity_v_Fhet.NCBI/"
            transrate_out_reverse = transratedir + sample + "_Fhet.NCBI_v_trinity/"
            transrate_assemblies_forward = transrate_out_forward + "assemblies.csv"
            transrate_assemblies_reverse = transrate_out_reverse + "assemblies.csv"
	    if os.path.isfile(transrate_assemblies_forward):
                data1 = parse_transrate_stats(transrate_assemblies_forward)
                data_frame1 = build_DataFrame(data_frame1, data1)
            else:  
                print("Running transrate forward...")
                transrate_forward(transratedir,transrate_out_forward,trinity_fasta,sample,reference)
            if os.path.isfile(transrate_assemblies_reverse):
                data2 = parse_transrate_stats(transrate_assemblies_reverse)
                data_frame2 = build_DataFrame(data_frame2,data2)
            else:
                print("Running transrate reverse...")
                transrate_reverse(transratedir,transrate_out_reverse,trinity_fasta,sample,reference)
    return data_frame1, data_frame2

assemblydir = "/home/ljcohen/public_html/killifish/transcriptome_assemblies/"
#transratedir = "/home/ljcohen/osmotic_transrate/kfish2rae5g_mrna.combined/"
transratedir = "/home/ljcohen/osmotic_transrate/NCBI/"
#reference = "/home/ljcohen/reference/kfish2rae5/kfish2rae5g.mrna.combined"
reference = "/home/ljcohen/reference/Fhet_ncbi/rna.fa"
clusterfunc_py3.check_dir(transratedir)
listoffiles = os.listdir(assemblydir)
data_frame1 = pd.DataFrame()
data_frame2 = pd.DataFrame()
data_frame1,data_frame2 = execute(data_frame1,data_frame2,listoffiles, assemblydir, transratedir,reference)
data_frame1.to_csv("../evaluation_data/transrate_reference_trinity_v_Fhet.NCBI.csv")
data_frame2.to_csv("../evaluation_data/transrate_reference_Fhet.NCBI_v_trinity.csv")
print("Forward transrate reference written: evaluation_data/transrate_reference_trinity_v_Fhet.NCBI.csv")
print("Reverse transrate reference written: evaluation_data/transrate_reference_Fhet.NCBI_v_trinity.csv")
