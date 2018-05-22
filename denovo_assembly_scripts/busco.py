import os
import os.path
from os.path import basename
# custom Lisa module
import clusterfunc_py3

def run_busco(fasta,sample,buscodir):
    busco_command="""
run_BUSCO.py -i {} -o {} -l ~/bin/busco/metazoa_odb9 -m tran --cpu 8
""".format(fasta,sample)
    print(busco_command)
    commands = [busco_command]
    process_name = "busco"
    module_name_list = ""
    filename = sample
    clusterfunc_py3.sbatch_file(buscodir, process_name,module_name_list, filename, commands)

def execute(assemblies,assemblydir,buscodir):
    for assembly in assemblies:
        if assembly.endswith(".fasta"):
            fasta = assemblydir + assembly
            sample = assembly.split(".")[0]
            run_busco(fasta,sample,buscodir)
 
buscodir = "/home/ljcohen/osmotic_busco_v3/metazoa/"
assemblydir = "/home/ljcohen/public_html/killifish/transcriptome_assemblies/"
listoffiles = os.listdir(assemblydir)
execute(listoffiles,assemblydir,buscodir)
