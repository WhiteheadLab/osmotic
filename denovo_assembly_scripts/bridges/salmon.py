import os
import os.path
import subprocess
from subprocess import Popen, PIPE
# custom Lisa module
import clusterfunc

def salmon_index(salmondir,assembly,genus_species):
    salmon_index_string="""
#source /home/ljcohen/.bashrc
cd {}
salmon index --index {} --transcripts {} --type quasi --threads 8
""".format(salmondir,genus_species,assembly)
    print(salmon_index_string)
    commands=[salmon_index_string]
    process_name="salmon_index_"+genus_species
    module_name_list=""
    #s=subprocess.Popen(salmon_index_string,shell=True)
    #s.wait()
    clusterfunc.sbatch_file(salmondir,process_name,module_name_list,genus_species,commands)
    #return salmon_index_string

def quant_salmon(salmondir,trimdir,genus_species,indexdir):
    salmon_string="""
source /home/ljcohen/.bashrc
cd {}
for i in {}{}*trim_1P.fq
do
        BASE=$(basename $i .trim_1P.fq)
        salmon quant -i {}{} -p 12 --validateMappings --libType IU -1 {}$BASE.trim_1P.fq -2 {}$BASE.trim_2P.fq -o $BASE.quant
done
""".format(salmondir,trimdir,genus_species,indexdir,genus_species,trimdir,trimdir)
    #salmon_index_string  = salmon_index(salmondir,genus_species)
    commands=[salmon_string]
    process_name="salmon"
    module_name_list=""
    clusterfunc.sbatch_file(salmondir,process_name,module_name_list,genus_species,commands)

def execute(assemblies,salmondir,trimdir,indexdir):
    for assembly in assemblies:
        if assembly.endswith(".fasta"):
            genus_species = assembly.split(".")[0]
            print(genus_species)
            #salmon_index(indexdir,assembly,genus_species)
            quant_salmon(salmondir,trimdir,genus_species,indexdir)

indexdir="/pylon5/bi5fpmp/ljcohen/kfish_salmon_ref/"
trimdir="/pylon5/bi5fpmp/ljcohen/kfish_trimmed/"
salmondir="/pylon5/bi5fpmp/ljcohen/kfish_salmon/"
assemblies=os.listdir(indexdir)
#assemblies=["F_heteroclitusMDPP","F_heteroclitusMDPL"]
execute(assemblies,salmondir,trimdir,indexdir)
