import os
import os.path
import clusterfunc

def run_sourmash(ref_dir,ref):
    sourmash="""
sourmash compute --scaled 10000 {}{} -o {}ensembl.sig -k 31
""".format(ref_dir,ref,ref_dir)
    sourmash_command=[sourmash]
    module_load_list=[]
    process_name="sourmashcompute"
    clusterfunc.sbatch_file(ref_dir,process_name,module_load_list,"ensembl",sourmash_command)

def execute(ref_dir,ref):
    run_sourmash(ref_dir,ref)

ref_dir = "/pylon5/bi5fpmp/ljcohen/Fhet_reference_genome/ensembl/"
ref = "Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.dna.toplevel.fa"
execute(ref_dir,ref)
