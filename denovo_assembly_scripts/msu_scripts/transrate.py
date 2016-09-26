import os
import os.path
from os.path import basename
import subprocess
from subprocess import Popen, PIPE
import clusterfunc



def get_assemblies(assemblydir):
	genus_species_dirs=os.listdir(assemblydir)
	for genus_species in genus_species_dirs:
		left=assemblydir+genus_species+"/"+genus_species+".left.fq"
		right=assemblydir+genus_species+"/"+genus_species+".right.fq"
		if os.path.isfile(left):
			print left
		else:
			print "there's a problem:",left
		if os.path.isfile(right):
			print right
		else:
			print "there's a problem:",right
		trinity_out_dir=assemblydir+genus_species+"/trinity_out/"
		trinity_fasta=trinity_out_dir+"Trinity.fasta"
		if os.path.isfile(trinity_fasta):
			fixed_trinity_fasta=fix_fasta(trinity_fasta,trinity_out_dir,genus_species)
			transrate_command=transrate(fixed_trinity_fasta,genus_species,left,right)
        		transrate_command=[transrate_command]
			module_load_list=["BLAST+/2.2.31"]
        		process_name="transrate"
        		clusterfunc.qsub_file(trinity_out_dir,process_name,module_load_list,genus_species,transrate_command)
		else:
			print "Assembly not completed:",genus_species

def transrate(fixed_trinity_fasta,genus_species,left,right):
	transrate_command="""
transrate --assembly {} --left {} --right {} --reference /mnt/home/ljcohen/reference/kfish2evg367mixx11pub2.mrna --threads 32 > {}.transrate.out
""".format(fixed_trinity_fasta,left,right,genus_species)	
	print transrate_command
	return transrate_command



def fix_fasta(trinity_fasta,trinity_dir,sample):
	trinity_out=trinity_dir+sample+".Trinity.fixed.fa"
	fix="""
sed 's_|_-_g' {} > {}
""".format(trinity_fasta,trinity_out)
	#s=subprocess.Popen(fix,shell=True)
	print fix
        #s.wait()
	return trinity_out 

assemblydir="/mnt/research/ged/lisa/osmotic_killifish/assemblies_msu/"
get_assemblies(assemblydir)
