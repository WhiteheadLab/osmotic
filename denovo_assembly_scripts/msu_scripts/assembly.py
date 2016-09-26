import os
import os.path
from os.path import basename
import subprocess
from subprocess import Popen, PIPE
import clusterfunc


def run_trinity(assemblydir):
        assemblydirs=os.listdir(assemblydir)
        for genus_species in assemblydirs:
                print genus_species
		genus_species_dir=assemblydir+genus_species+"/"
                listoffiles=os.listdir(genus_species_dir)
		print listoffiles
		trinity_command="""
set -x
# stops execution if there is an error
set -e
if [ -f {}trinity_out/Trinity.fasta ]; then exit 0 ; fi
if [ -d {}trinity_out ]; then mv {}trinity_out_dir {}trinity_out_dir0 || true ; fi

Trinity --left {}{}.left.fq \\
--right {}{}.right.fq \\
--output {}trinity_out --seqType fq --monitoring --bflyCalculateCPU --max_memory 150G	\\
--CPU 32
""".format(genus_species_dir,genus_species_dir,genus_species_dir,genus_species_dir,genus_species_dir,genus_species,genus_species_dir,genus_species,genus_species_dir)
		print trinity_command
		trinity_command=[trinity_command]
		#module_load_list=["trinity/2.0.5"]
		module_load_list=["trinity/2.2.0"]
		process_name="trinity"
		clusterfunc.qsub_file(genus_species_dir,process_name,module_load_list,genus_species,trinity_command)

assemblydir="/mnt/research/ged/lisa/osmotic_killifish/fix_assemblies_msu/"
#clusterfunc.check_dir(assemblydir)
run_trinity(assemblydir)
