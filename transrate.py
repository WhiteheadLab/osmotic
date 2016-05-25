import os
import os.path
from os.path import basename
import subprocess
from subprocess import Popen, PIPE
import clusterfunc



def get_assemblies(assemblydir,readsdir):
	#genus_species_dirs=os.listdir(readsdir)
	genus_species_dirs=["F_diaphanus/F_diaphanus.trinity","F_sciadicus/F_sciadicus.trinity"]
	for dirname in genus_species_dirs:
		genus_species=dirname.split("/")[0]
		genus_species_dir=assemblydir+genus_species+"/"
		reads=readsdir+genus_species+"/"+"diginorm/"
		left=reads+genus_species+".left.fq"
		right=reads+genus_species+".right.fq"
		if os.path.isfile(left):
			print left
		else:
			print "there's a problem:",left
		if os.path.isfile(right):
			print right
		else:
			print "there's a problem:",right
		trinity_out_dir1=assemblydir+dirname+".1/"
		trinity_out_dir2=assemblydir+dirname+".2/"
		trinity_fasta1=trinity_out_dir1+"Trinity.fasta"
		trinity_fasta2=trinity_out_dir2+"Trinity.fasta"
		fixed_trinity_fasta1=fix_fasta(trinity_fasta1,trinity_out_dir1,genus_species)
		print fixed_trinity_fasta1
		fixed_trinity_fasta2=fix_fasta(trinity_fasta2,trinity_out_dir2,genus_species)
		print fixed_trinity_fasta2
		transrate_command=transrate(trinity_out_dir1,trinity_out_dir2,fixed_trinity_fasta1,fixed_trinity_fasta2,genus_species,left,right,assemblydir)
        	print transrate_command
		transrate_command=[transrate_command]
		module_load_list=["blast/2.2.29"]
        	process_name="transrate"
        	clusterfunc.sbatch_file(genus_species_dir,process_name,module_load_list,genus_species,transrate_command)

def transrate(trinity_out_dir1,trinity_out_dir2,fixed_trinity_fasta1,fixed_trinity_fasta2,genus_species,left,right,assemblydir):
	transrate_command="""
transrate --assembly {},{} \\
--left {} \\
--right {} \\
--reference /home/ljcohen/reference/kfish2rae5/kfish2rae5g.mrna.combined \\
--threads 32
""".format(fixed_trinity_fasta1,fixed_trinity_fasta2,left,right)	

#	transrate_command="""
#transrate --assembly {} \\
#--reference /home/ljcohen/reference/kf2evg367mixx11/kfish2evg367mixx11pub2.mrna \\
#--output {}kf2evg367mixx11.transrate/ \\
#--threads 32
#""".format(fixed_trinity_fasta,trinity_out_dir)
	return transrate_command

#--reference /home/ljcohen/reference/kfish2rae5/kfish2rae5g.mrna.combined
#--reference /home/ljcohen/reference/kf2evg367mixx11/kfish2evg367mixx11pub2.mrna \\

def fix_fasta(trinity_fasta,trinity_dir,sample):
	trinity_out=trinity_dir+sample+".Trinity.fixed.fa"
	fix="""
sed 's_|_-_g' {} > {}
""".format(trinity_fasta,trinity_out)
	#s=subprocess.Popen(fix,shell=True)
	print fix
        #s.wait()
	return trinity_out 

assemblydir="/home/ljcohen/msu_assemblies_finished/"
readsdir="/home/ljcohen/osmotic_trim_"
get_assemblies(assemblydir,readsdir)
