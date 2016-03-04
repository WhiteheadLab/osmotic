import os
import os.path
from os.path import basename
import subprocess
from subprocess import Popen, PIPE
import clusterfunc


def run_trinity(assemblydir):
        assemblydirs=os.listdir(assemblydir)
        for genus_species in assemblydirs:
                genus_species_dir=assemblydir+genus_species+"/"
                listoffiles=os.listdir(genus_species_dir)
		trinity_command="""
set -x
# stops execution if there is an error
set -e
if [ -f {}trinity_out/Trinity.fasta ]; then exit 0 ; fi
if [ -d {}trinity_out ]; then mv {}trinity_out_dir {}trinity_out_dir0 || true ; fi

Trinity --left {}{}.left.fq \\
--right {}{}.right.fq \\
--output {}trinity_out --seqType fq --max_memory 14G	\\
--CPU ${{THREADS:-2}}
""".format(genus_species_dir,genus_species_dir,genus_species_dir,genus_species_dir,genus_species_dir,genus_species,genus_species_dir,genus_species,genus_species_dir)
		print trinity_command
		trinity_command=[trinity_command]
		module_load_list=["rsem/1.2.23","trinity/2.0.5"]
		process_name="trinity"
		clusterfunc.sbatch_file(genus_species_dir,process_name,module_load_list,genus_species,trinity_command)

assemblydir="/home/ljcohen/osmotic_trim/assemblies/"
clusterfunc.check_dir(assemblydir)
run_trinity(assemblydir)
