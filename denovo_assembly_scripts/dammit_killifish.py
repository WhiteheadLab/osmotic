import os
import os.path
from os.path import basename
import subprocess
from subprocess import Popen, PIPE
import clusterfunc_py3


def dammit(trinity_fasta,genus_species,dammitdir):
    dammit_command="""
source /home/ljcohen/.bashrc
source activate /home/ljcohen/anaconda2/envs/py3.dammit

mkdir /scratch/$SLURM_JOBID
cd /scratch/$SLURM_JOBID

dammit annotate {} --busco-group metazoa --user-databases /home/ljcohen/reference/kfish2rae5/kfish2rae5g.combined.aa -o {}.dammit_out --n_threads 12

cp -r /scratch/$SLURM_JOBID/{}.dammit_out {}
rm -rf /scratch/$SLURM_JOBID*
""".format(trinity_fasta,genus_species,genus_species,dammitdir)
    return dammit_command

def run_dammit(dammit_command,dammitdir,genus_species):
    dammitstring=[dammit_command]
    process_name="dammit"
    module_name_list=""
    clusterfunc_py3.sbatch_file(dammitdir,process_name,module_name_list,genus_species,dammitstring)

def execute(assemblies,assemblydir,dammitdir):
    for assembly in assemblies:
        if assembly.endswith(".trinity_out.Trinity.fasta"):
            genus_species = assembly.split(".")[0]
            print(genus_species)
            trinity_fasta=assemblydir+assembly
            dammit_command=dammit(trinity_fasta,genus_species,dammitdir)
            run_dammit(dammit_command,dammitdir,genus_species)

assemblydir="/home/ljcohen/osmotic_assemblies_farm/"
dammitdir = "/home/ljcohen/osmotic_dammit_kfish2rae5/"
clusterfunc_py3.check_dir(dammitdir)
assemblies=os.listdir(assemblydir)
execute(assemblies,assemblydir,dammitdir)
