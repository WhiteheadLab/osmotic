import os
import os.path
from os.path import basename
# custom Lisa module
import clusterfunc_py3

def run_busco(fasta,sample,buscodir):
    busco_command="""
module unload python
source /home/ljcohen/.bashrc
conda activate dammit_new
cd /mnt/scratch/ljcohen/osmotic_killifish/busco
run_BUSCO.py -i {} -o {}.actino -l ~/reference/actinopterygii_odb9 -m tran --cpu 16
""".format(fasta,sample)
    print(busco_command)
    commands = [busco_command]
    process_name = "busco"
    module_name_list = ""
    filename = sample
    clusterfunc_py3.qsub_file(buscodir, process_name,module_name_list, filename, commands)

def execute(assemblies,assemblydir,buscodir):
    for assembly in assemblies:
        if assembly.endswith(".fasta"):
            fasta = assemblydir + assembly
            sample = assembly.split(".")[0]
            run_busco(fasta,sample,buscodir)
 
buscodir = "/mnt/scratch/ljcohen/osmotic_killifish/busco/"
assemblydir = "/home/ljcohen/osmotic_killifish_assemblies/assemblies/"
listoffiles = os.listdir(assemblydir)
execute(listoffiles,assemblydir,buscodir)
