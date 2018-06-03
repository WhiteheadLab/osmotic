import os
import os.path
from os.path import basename
# custom Lisa module
import clusterfunc_py3
import pandas as pd

def get_pairs(listoffiles,basedir):
    pairs={}
    for reads_filename in listoffiles:
        if reads_filename.endswith(".fq"):
            print(reads_filename)
            filename=basedir+reads_filename
            fields=reads_filename.split("_")
            genus = fields[0]
            species = fields[1]
            genus_species = genus + "_" + species
            #print(pairs)
            if genus_species in pairs:
                pairs[genus_species].append(filename)
                sort_pairs = pairs[genus_species]
                sort_pairs.sort()
                pairs[genus_species] = sort_pairs
            else:
                pairs[genus_species]=[filename]
    return pairs

def combined_trim_treads(pairs):
    for sample in pairs:
       files = " ".join(pairs[sample])
       combine_string = "cat "+files+" > "+sample+".fq"
       print(combine_string)
    
def fix_fasta(trinity_fasta, fixed_trinity_fasta):
    fix = """
sed 's_|_-_g' {} > {}
""".format(trinity_fasta, fixed_trinity_fasta)
    #s=subprocess.Popen(fix,shell=True)
    print(fix)
    #s.wait()

def transrate(transratedir,transrate_out,trinity_fasta,sample,left,right):
    transrate_command = """
source ~/.bashrc
module load GNU/4.8.3
module unload python
module load parallel
source activate dammit_new
transrate --assembly={} --threads=16 \
--left={} \
--right={} \
--output={}
""".format(trinity_fasta,left,right,transrate_out)
    print(transrate_command)
    commands = [transrate_command]
    process_name = "transrate"
    module_name_list = ""
    filename = sample
    clusterfunc_py3.qsub_file(transratedir, process_name,module_name_list, filename, commands)

def parse_transrate_stats(transrate_assemblies):
    print(transrate_assemblies)
    if os.stat(transrate_assemblies).st_size != 0:
        data = pd.DataFrame.from_csv(transrate_assemblies, header=0, sep=',')
        return data

def build_DataFrame(data_frame, transrate_data):
    # columns=["n_bases","gc","gc_skew","mean_orf_percent"]
    frames = [data_frame, transrate_data]
    data_frame = pd.concat(frames)
    return data_frame

def execute(data_frame, listoffiles, listofassemblies, assemblydir, transratedir):
    pairs_dictionary=get_pairs(listoffiles,basedir)
    # construct an empty pandas dataframe to add on each assembly.csv to
    for assembly in listofassemblies:
        print(assembly)
        if assembly.endswith(".fasta"):
            sample = assembly.split("_")[0] + "_" + assembly.split("_")[1][:-8]
            if sample in pairs_dictionary:
                files = pairs_dictionary[sample]
                left = files[0]
                right = files[1]
                trinity_fasta = assemblydir + assembly
                transrate_out = transratedir + sample + "/"
                transrate_assemblies = transrate_out + "/" + "assemblies.csv"
                if os.path.isfile(transrate_assemblies):
                    data = parse_transrate_stats(transrate_assemblies)
                    data_frame = build_DataFrame(data_frame, data)
                else: 
                    print("Running transrate...")
                    #transrate(transratedir,transrate_out,trinity_fasta,sample,left,right)
    return data_frame

assemblies = "/mnt/home/ljcohen/osmotic_killifish_assemblies/assemblies/"
basedir = "/mnt/scratch/ljcohen/osmotic_killifish/combined_trimmed/"
transratedir = "/mnt/scratch/ljcohen/osmotic_killifish/transrate/"
clusterfunc_py3.check_dir(transratedir)
listoffiles = os.listdir(basedir)
listofassemblies = os.listdir(assemblies)
data_frame = pd.DataFrame()
data_frame = execute(data_frame, listoffiles, listofassemblies, assemblies, transratedir)
data_frame.to_csv("../evaluation_data/transrate_scores.csv")
