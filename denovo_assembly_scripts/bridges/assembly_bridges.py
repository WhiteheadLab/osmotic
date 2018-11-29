import os
import os.path
from os.path import basename
import subprocess
from subprocess import Popen, PIPE
import clusterfunc


def run_trinity(assemblydir,filesdir,species):
    trinity_command="""
# copy files over to the $LOCAL node storage, associated with the running job
# (because Trinity does a lot of read/writes (IO-bound) 

SPECIES={}
PROJECTDIR=$LOCAL/$SPECIES/
mkdir $PROJECTDIR
cd $PROJECTDIR
cp {}{}.left.fq .
cp {}{}.right.fq .

# run Trinity

Trinity --left {}.left.fq \\
--right {}.right.fq \\
--output $PROJECTDIR/$SPECIES.trinity_out \\
--full_cleanup --seqType fq \\
--max_memory 1500G \\
--CPU 42

# Grab the assembly file and copy it to your storage

cp $PROJECTDIR/$SPECIES.trinity_out.Trinity.fasta {}

# remove temp files on $LOCAL
rm -rf $PROJECTDIR
""".format(species,filesdir,species,filesdir,species,species,species,assemblydir)
    trinity_command=[trinity_command]
    module_load_list=["samtools/1.7","jellyfish2/2.2.6","bowtie2/2.3.4.1","salmon/0.9.1","trinity/2.8.4"]
    process_name="trinity"
    clusterfunc.sbatch_file(assemblydir,process_name,module_load_list,species,trinity_command)

def execute(assemblydir,filesdir,txome_files):
    species_repeat=[]
    for filename in txome_files:
        if filename.endswith(".fq"):
            species=filename.split(".")[0]
            if species not in species_repeat:
                run_trinity(assemblydir,filesdir,species)
                species_repeat.append(species)
assemblydir="/pylon5/bi5fpmp/ljcohen/kfish_trinity/"
filesdir="/pylon5/bi5fpmp/ljcohen/kfish_txomes/"
clusterfunc.check_dir(assemblydir)
txome_files = os.listdir(filesdir)
execute(assemblydir,filesdir,txome_files)
