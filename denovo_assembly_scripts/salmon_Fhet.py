import os
import os.path
import subprocess
from subprocess import Popen, PIPE
# custom Lisa module
import clusterfunc

def salmon_index(salmondir,genus_species,trinity_fasta):
	index=genus_species+"_index"
	salmon_index_string="""
cd {}
salmon index --transcripts /home/ljcohen/reference/kfish2rae5/kfish2rae5g.main.pub.mrna --index kfish2rae5g.main.tub.mrna_index --type quasi --threads 8
""".format(salmondir,index,trinity_fasta)
	return salmon_index_string,index

def quant_salmon(salmondir,trimdir,genus_species,trinity_fasta):
	salmon_string="""
cd {}
for i in {}{}*trim_1P.fq
do
	BASE=$(basename $i .trim_1P.fq)
        salmon quant -i kfish2rae5g.main.tub.mrna_index --libType IU -1 {}$BASE.trim_1P.fq -2 {}$BASE.trim_2P.fq -o $BASE.quant
done
""".format(salmondir,trimdir,genus_species,trimdir,trimdir)
	salmonstring=[salmon_string]
        process_name="salmon"
        module_name_list=""
        clusterfunc.sbatch_file(salmondir,process_name,module_name_list,genus_species,salmonstring)

def execute(assemblies,salmondir,assemblydir,basedir,trimdir):
	for assembly in assemblies:
		genus_species = assembly.split(".")[0]
		print genus_species
		trinity_fasta=assemblydir+assembly
		quant_salmon(salmondir,trimdir,genus_species,trinity_fasta)


basedir="/home/ljcohen/osmotic/"
trimdir="/home/ljcohen/osmotic_trim/"
assemblydir="/home/ljcohen/osmotic_assemblies_farm/"
salmondir="/home/ljcohen/osmotic_salmon_Fhet_cds/"
clusterfunc.check_dir(salmondir)
assemblies=os.listdir(assemblydir)
execute(assemblies,salmondir,assemblydir,basedir,trimdir)
