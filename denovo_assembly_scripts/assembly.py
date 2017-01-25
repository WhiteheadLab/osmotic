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
PROJECTDIR={}
SPECIES={}

mkdir /scratch/$SLURM_JOBID
cd /scratch/$SLURM_JOBID
cp $PROJECTDIR/*.fq .
# monitor resources used
sar -u -r -d -o /scratch/$SLURM_JOBID/$SPECIES.times.dat 1

Trinity --left $SPECIES.left.fq \\
--right $SPECIES.right.fq \\
--output /scratch/$SLURM_JOBID/$SPECIES.trinity_out \\
--full_cleanup --seqType fq \\
--max_memory 256G \\
--CPU 30
cp /scratch/$SLURM_JOBID/$SPECIES.trinity_out.Trinity.fasta /home/ljcohen/osmotic_assemblies_farm/
# extract data from resource monitoring
sar -d -p -f /scratch/$SLURM_JOBID/$SPECIES.times.dat > /home/ljcohen/osmotic_assemblies_farm/$SPECIES.disk.txt
sar -u -f /scratch/$SLURM_JOBID/$SPECIES.times.dat > /home/ljcohen/osmotic_assemblies_farm/$SPECIES.cpu.txt
sar -r -f /scratch/$SLURM_JOBID/$SPECIES.times.dat > /home/ljcohen/osmotic_assemblies_farm/$SPECIES.ram.txt
gzip /home/ljcohen/osmotic_assemblies_farm/$SPECIES.*.txt
# remove temp files on scratch
rm -rf /scratch/$SLURM_JOBID*
""".format(genus_species_dir,genus_species)
		print trinity_command
		trinity_command=[trinity_command]
		module_load_list=["rsem/1.2.23","trinity/2.2.0","java/1.8"]
		process_name="trinity"
		clusterfunc.sbatch_file(genus_species_dir,process_name,module_load_list,genus_species,trinity_command)

assemblydir="/home/ljcohen/assemblies/"
clusterfunc.check_dir(assemblydir)
run_trinity(assemblydir)
