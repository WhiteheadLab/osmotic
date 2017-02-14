import os
import os.path
from os.path import basename
import subprocess
from subprocess import Popen, PIPE
import clusterfunc_py3
import pandas as pd

def get_assemblies(assemblydir,transrate_outdir1,transrate_outdir2):
    assemblyfiles=os.listdir(assemblydir)
    #genus_species_dirs=["F_heteroclitus.MDPL","F_heteroclitus.MDPP"]
    for assembly in assemblyfiles:
        if assembly.endswith(".dammit.fasta"):
            genus_species=assembly.split(".")[0]
            trinity_fasta=assemblydir + assembly
            fixed_trinity_fasta=fix_fasta(trinity_fasta,assemblydir,genus_species)
            print(fixed_trinity_fasta)
            transrate_reference_forward(transrate_outdir1,fixed_trinity_fasta,genus_species)
            transrate_reference_forward(transrate_outdir2,fixed_trinity_fasta,genus_species)

def transrate_reference_forward(transrate_outdir1,fixed_trinity_fasta,genus_species):
    transrate_command="""
transrate --assembly {} \\
--reference /home/ljcohen/reference/kfish2rae5/kfish2rae5g.main.pub.mrna \\
--output {}{}_v_Fhet \\
--threads 16
""".format(fixed_trinity_fasta,transrate_outdir1,genus_species)
    print(transrate_command)
    transrate_command=[transrate_command]
    module_load_list=["blast/2.2.29"]
    process_name="transrate"
    #clusterfunc_py3.sbatch_file(transrate_outdir1,process_name,module_load_list,genus_species,transrate_command)

def transrate_reference_reverse(transrate_outdir2,fixed_trinity_fasta,genus_species):
    transrate_command="""
transrate --assembly=/home/ljcohen/reference/kf2evg367mixx11/kfish2evg367mixx11pub2.mrna \\
--reference={} \\
--output {}Fhet_v_{} \\
--threads 16
""".format(fixed_trinity_fasta,transrate_outdir2,genus_species)
    print(transrate_command)
    transrate_command=[transrate_command]
    module_load_list=["blast/2.2.29"]
    process_name="transrate_rev"
    #clusterfunc_py3.sbatch_file(transrate_outdir2,process_name,module_load_list,genus_species,transrate_command)

def parse_transrate_stats(transrate_assemblies):
    data = pd.DataFrame.from_csv(transrate_assemblies, header=0, sep=',')
    return data

def build_DataFrame(data_frame, transrate_data):
    # columns=["n_bases","gc","gc_skew","mean_orf_percent"]
    frames = [data_frame, transrate_data]
    data_frame = pd.concat(frames)
    return data_frame

def get_contigs_data(data_frame1, data_frame2,transrate_dir1,transrate_dir2):
    listofdirs1 = os.listdir(transrate_dir1)
    listofdirs2 = os.listdir(transrate_dir2)
    for dirname1 in listofdirs1:
        if dirname1 != "sbatch_files":
            transrate_dirname1 = transrate_dir1 + dirname1 + "/"
            transrate_assemblies1 = transrate_dirname1 + "assemblies.csv"
            if os.path.isfile(transrate_assemblies1):
                print(transrate_assemblies1)
                data1 = parse_transrate_stats(transrate_assemblies1)
                data_frame1 = build_DataFrame(data_frame1, data1)
            else:
                print("File missing:", transrate_assemblies1)
    for dirname2 in listofdirs2:
        if dirname2 != "sbatch_files":
            transrate_dirname2 = transrate_dir2 + dirname2 + "/"
            transrate_assemblies2 = transrate_dirname2 + "assemblies.csv"
            if os.path.isfile(transrate_assemblies2):
                print(transrate_assemblies2)
                data2 = parse_transrate_stats(transrate_assemblies2)
                data_frame2 = build_DataFrame(data_frame2,data2)
            else:
                print("Files missing:",transrate_assemblies2)
    return data_frame1,data_frame2

#--reference /home/ljcohen/reference/kfish2rae5/kfish2rae5g.mrna.combined
#--reference /home/ljcohen/reference/kf2evg367mixx11/kfish2evg367mixx11pub2.mrna \\

def fix_fasta(trinity_fasta,trinity_dir,sample):
    trinity_out=trinity_dir+sample+".fixed.fasta"
    fix="""
sed 's_|_-_g' {} > {}
""".format(trinity_fasta,trinity_out)
    #s=subprocess.Popen(fix,shell=True)
    print(fix)
    #s.wait()
    return trinity_out 

assemblydir="/home/ljcohen/osmotic_assemblies_farm/"
transrate_outdir1 = "/home/ljcohen/transrate_reference_osmotic_v_Fhet/"
transrate_outdir2 = "/home/ljcohen/transrate_reference_Fhet_v_osmotic/"
#get_assemblies(assemblydir,transrate_outdir1,transrate_outdir2)
data_frame1 = pd.DataFrame()
data_frame2 = pd.DataFrame()
data_frame1, data_frame2 = get_contigs_data(data_frame1, data_frame2,transrate_outdir1,transrate_outdir2)
data_frame1.to_csv("../evaluation_data/transrate_reference_trinity_v_Fhet.csv")
data_frame2.to_csv("../evaluation_data/transrate_reverse_Fhet_v_trinity.csv")
print("Reference scores written.")
print("Reverse reference scores written.")
