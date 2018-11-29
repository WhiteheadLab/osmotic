import os
import os.path
import clusterfunc

def run_sourmash(reads_dir,read_file,species,sourmash_dir):
    sourmash="""
sourmash compute --scaled 10000 {}{} -o {}{}.sig -k 31
""".format(reads_dir,read_file,sourmash_dir,species)
    sourmash_command=[sourmash]
    module_load_list=[]
    process_name="sourmashcompute"
    clusterfunc.sbatch_file(reads_dir,process_name,module_load_list,species,sourmash_command)

def execute(reads_dir,sourmash_dir):
    species_repeat = []
    reads_files = os.listdir(reads_dir)
    for filename in reads_files:
        if filename.endswith("left.fq"):
            species=filename.split(".")[0]
            if species not in species_repeat:
                run_sourmash(reads_dir,filename,species,sourmash_dir)
                species_repeat.append(species)

reads_dir = "/pylon5/bi5fpmp/ljcohen/kfish_txomes/"
sourmash_dir = "/pylon5/bi5fpmp/ljcohen/kfish_sourmash_reads/"
clusterfunc.check_dir(sourmash_dir)
execute(reads_dir,sourmash_dir)
