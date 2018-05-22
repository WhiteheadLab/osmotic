import os
# custom Lisa module
import clusterfunc_msu


def get_dammit_string(trinity_fasta):
    j = """
source /mnt/home/ljcohen/.bashrc
module load GNU/4.8.3
module unload python
module load parallel
source activate dammit_new
module load LAST/737
export DAMMIT_DB_DIR=/mnt/research/ged/data/dammitdb/
dammit annotate {} --busco-group eukaryota --user-databases /mnt/home/ljcohen/reference/kfish2rae5g.aa.fasta --n_threads 14
""".format(trinity_fasta)
    return j

def run_dammit(dammit_string,dammitdir,mmetsp):
    dammit_command = [dammit_string]
    process_name = "dammit"
    module_name_list = []
    filename = mmetsp
    clusterfunc_msu.qsub_file(dammit_dir, process_name, module_name_list, filename, dammit_command)

def execute(assemblies, basedir, dammit_dir):
    for assembly in assemblies:
        if assembly.endswith(".fasta"):
            mmetsp = assembly.split(".")[0]
            print(mmetsp)
            trinity_fasta = basedir+assembly
            print(trinity_fasta)
            dammit_string = get_dammit_string(trinity_fasta)
            run_dammit(dammit_string, dammit_dir, mmetsp)

basedir = "/mnt/home/ljcohen/killifish_assemblies/"
dammit_dir = "/mnt/home/ljcohen/killifish_assemblies/"
assemblies = os.listdir(basedir)
execute(assemblies, basedir,dammit_dir)
