import os
import os.path
from os.path import basename
import subprocess
from subprocess import Popen, PIPE
import clusterfunc_py3


def dammit(trinity_fasta,genus_species,dammitdir):
    dammit_command="""
dammit annotate {} --busco-group metazoa --user-databases /home/ljcohen/osmotic/kfish_reference_genome_files/protein.fa -o {}{}.dammit_out --n_threads 16
""".format(trinity_fasta,dammitdir,genus_species)
    return dammit_command

def run_dammit(dammit_command,dammitdir,genus_species):
    dammitstring=[dammit_command]
    process_name="dammit_"+genus_species
    module_name_list=""
    clusterfunc_py3.sbatch_file(dammitdir,process_name,module_name_list,genus_species,dammitstring)

def execute(assemblies,assemblydir,dammitdir):
    for assembly in assemblies:
        genus_species = assembly.split(".")[0]
        print(genus_species)
        trinity_fasta=assemblydir+assembly
        dammit_command=dammit(trinity_fasta,genus_species,dammitdir)
        run_dammit(dammit_command,dammitdir,genus_species)

assemblydir="/home/ljcohen/osmotic_assemblies_farm/"
dammitdir = "/home/ljcohen/osmotic_damit/"
clusterfunc_py3.check_dir(dammitdir)
assemblies=os.listdir(assemblydir)
execute(assemblies,assemblydir,dammitdir)

