import os
import os.path
import subprocess
from subprocess import Popen, PIPE
# custom Lisa module
import clusterfunc_py3

def salmon_index(salmondir,genus_species):
	salmon_index_string="""
cd {}
salmon index --transcripts --index {} --type quasi --threads 8
""".format(salmondir,genus_species)
	return salmon_index_string

def quant_salmon(salmondir,trimdir,genus_species,indexdir):
        salmon_string="""
source /home/ljcohen/.bashrc
cd {}
for i in {}{}*trim_1P.fq
do
        BASE=$(basename $i .trim_1P.fq)
        salmon quant -i {}kfish2rae5g --libType IU -1 {}$BASE.trim_1P.fq -2 {}$BASE.trim_2P.fq -o $BASE.quant
done
""".format(salmondir,trimdir,genus_species,indexdir,trimdir,trimdir)
	#salmon_index_string  = salmon_index(salmondir,genus_species)
        commands=[salmon_string]
        process_name="salmon"
        module_name_list=""
        clusterfunc_py3.sbatch_file(salmondir,process_name,module_name_list,genus_species,commands)

def execute(assemblies,salmondir,assemblydir,trimdir,indexdir):
        for assembly in assemblies:
                genus_species = assembly.split(".")[0]
                print(genus_species)
                quant_salmon(salmondir,trimdir,genus_species,indexdir)

indexdir="/home/ljcohen/reference/kfish2rae5/"
trimdir="/home/ljcohen/osmotic_trim/"
assemblydir="/home/ljcohen/public_html/killifish/transcriptome_assemblies/"
salmondir="/home/ljcohen/osmotic_salmon_Fhet_reference/"
clusterfunc_py3.check_dir(salmondir)
assemblies=os.listdir(assemblydir)
#assemblies=["F_heteroclitusMDPP","F_heteroclitusMDPL"]
execute(assemblies,salmondir,assemblydir,trimdir,indexdir)
