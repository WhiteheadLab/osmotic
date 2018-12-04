import os
import os.path
from os.path import basename
import subprocess
from subprocess import Popen, PIPE
import clusterfunc


def run_plass(assemblydir,filesdir,species):
    plass_command="""
source activate plass
plass assemble {}{}.left.fq {}{}.right.fq {}{}_plass.fasta {}_tmp
""".format(filesdir,species,filesdir,species,assemblydir,species,species)
    plass_command=[plass_command]
    module_load_list=[]
    process_name="plass"
    clusterfunc.sbatch_file(assemblydir,process_name,module_load_list,species,plass_command)

def execute(assemblydir,filesdir,txome_files):
    species_repeat=[]
    for filename in txome_files:
        if filename.endswith(".fq"):
            species=filename.split(".")[0]
            if species not in species_repeat:
                run_plass(assemblydir,filesdir,species)
                species_repeat.append(species)
assemblydir="/pylon5/bi5fpmp/ljcohen/kfish_plass/"
filesdir="/pylon5/bi5fpmp/ljcohen/kfish_txomes/"
clusterfunc.check_dir(assemblydir)
txome_files = os.listdir(filesdir)
execute(assemblydir,filesdir,txome_files)
