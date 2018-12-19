import os
import os.path
from os.path import basename
import subprocess
from subprocess import Popen, PIPE
import clusterfunc

def run_dammit(aa_databases,assembly,species,databases,basedir,dammit_dir):
    dammit_string="""
source /home/ljcohen/.bashrc
source activate dammit_conda
export DAMMIT_DB_DIR=/pylon5/bi5fpmp/ljcohen/dammit
SPECIES={}
PROJECTDIR=$LOCAL/$SPECIES
mkdir $PROJECTDIR
cd $PROJECTDIR
cp {}{} .
cp -r {} .
dammit annotate {} --busco-group eukaryota \\
     --user-databases {} \\
     --output-dir {}.dammit \\
     --n_threads 14
cp -r {}.dammit {}
cd ..
rm -rf $PROJECTDIR
""".format(species,basedir,assembly,aa_databases,assembly,databases,species,species,dammit_dir)
    dammit_command = [dammit_string]
    process_name = "dammit" + "_" + species
    module_name_list = []
    filename = species
    clusterfunc.sbatch_file(dammit_dir, process_name, module_name_list, filename, dammit_command)

def execute(fasta_files,databases,basedir,dammit_dir,aa_databases):
    for assembly in fasta_files:
        if assembly.endswith(".fasta"):
            species=assembly.split(".")[0]
            print(species)
            run_dammit(aa_databases,assembly,species,databases,basedir,dammit_dir)

#dammit_dir="/pylon5/bi5fpmp/ljcohen/kfish_dammit/"
#dammit_dir="/pylon5/bi5fpmp/ljcohen/kfish_dammit_evigene/"
dammit_dir = "/pylon5/bi5fpmp/ljcohen/kfish_dammit_ncbi"
basedir = "/pylon5/bi5fpmp/ljcohen/kfish_trinity/"
aa_databases = "/pylon5/bi5fpmp/ljcohen/Fhet_reference_genome/"
ncbi = "Fhet_reference_genome/ncbi/protein.fa"
evigene = "Fhet_reference_genome/evigene/kfish2rae5g.pub.aa"
ensembl = "Fhet_reference_genome/ensembl/Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.pep.all.fa"
#databases = ncbi + " " + evigene + " " + ensembl
#databases = ensembl
#databases = evigene
databases = ncbi
fasta_files = os.listdir(basedir)
execute(fasta_files,databases,basedir,dammit_dir,aa_databases) 

