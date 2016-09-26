import os
import os.path
from os.path import basename
import subprocess
from subprocess import Popen, PIPE
import clusterfunc

def get_assemblies(assemblydir):
	genus_species_dirs=os.listdir(assemblydir)
	for genus_species in genus_species_dirs:
		trinity_out_dir=assemblydir+genus_species+"/"
		trinity_fasta=trinity_out_dir+genus_species+".Trinity.fixed.fa"
		if os.path.isfile(trinity_fasta):
			dammit_command=[dammit(trinity_fasta)]
			module_load_list=["blast/2.2.29"]
        		process_name="dammit"
        		clusterfunc.sbatch_file(trinity_out_dir,process_name,module_load_list,genus_species,dammit_command)
		else:
			print "Assembly not completed:",genus_species

def dammit(trinity_fasta):
	dammit_command="""
dammit annotate {} --busco-group eukaryota --database-dir /home/ljcohen/reference/dammit --n_threads 16
""".format(trinity_fasta)
	return dammit_command


assemblydir="/home/ljcohen/msu_assemblies_finished/"
get_assemblies(assemblydir)
