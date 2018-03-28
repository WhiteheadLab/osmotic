import os
import os.path
import clusterfunc_py3

def run_trim_low_abund(trim_filename,trim_low_abund_dir):
#    trim_low_abund="""
#trim-low-abund.py -k 21 -C 2 --variable-coverage {}
#""".format(trim_filename)
    trim_low_abund="""
trim-low-abund.py -k 21 -C 2 {}
""".format(trim_filename)
    commands = [trim_low_abund]
    process_name = "trim-low-abund"
    module_name_list = ""
    filename = genus_species_sample
    clusterfunc_py3.sbatch_file(trim_low_abund_dir,process_name, module_name_list, filename, commands)

def head_reads(sourmash_filename,sourmash_dir,genus_species_sample):
    head_command="""
head -4000000 {} > /home/ljcohen/osmotic_sourmash/{}.head
""".format(sourmash_filename,genus_species_sample)
    commands = [head_command]
    process_name = "head"
    module_name_list = ""
    filename = genus_species_sample
    clusterfunc_py3.sbatch_file(sourmash_dir,process_name, module_name_list, filename, commands)

def get_sourmash_command(sourmash_filename,sourmash_dir,genus_species_sample):
#    sourmash_command="""
#sourmash compute --force -k 21 --track-abundance --scaled 20000 {}
#""".format(sourmash_filename)
    sourmash_command="""
sourmash compute --force -k 21 --scaled 20000 {}
""".format(sourmash_filename)
    commands = [sourmash_command]
    process_name = "sourmash"
    module_name_list = ""
    filename = genus_species_sample
    print(sourmash_command)
    clusterfunc_py3.sbatch_file(sourmash_dir,process_name, module_name_list, filename, commands)

def get_sourmash_streaming(sourmash_filename,sourmash_dir,genus_species_sample):
    sourmash_command="""
head  -n 4000000 {} | sourmash compute -k 51 --name {} --scaled 500 --dna - -o {}{}.sig
""".format(sourmash_filename,genus_species_sample[:-10],sourmash_dir,genus_species_sample)
    commands = [sourmash_command]
    process_name = "sourmash"
    module_name_list = ""
    filename = genus_species_sample
    print(sourmash_command)
    clusterfunc_py3.sbatch_file(sourmash_dir,process_name, module_name_list, filename, commands)

#trim_dir = "/home/ljcohen/osmotic_trim/"
trim_low_abund_dir = "/home/ljcohen/osmotic_trim/trim_low_abund/"
#trim_dir = "/home/ljcohen/mislabeled_files_removed_osmotic/"
#trim_dir = "/home/ljcohen/osmotic_trim/trim_low_abund/sbatch_files/"
trim_dir = "/home/ljcohen/osmotic_trim/trim_low_abund/"
#trim_dir = "/home/ljcohen/osmotic_sourmash/"
sourmash_dir = "/home/ljcohen/osmotic_sourmash/"
trim_files = os.listdir(trim_dir)
completed_sourmash = "/home/ljcohen/osmotic_sourmash/sbatch_files/"
complete_sourmash_files = os.listdir(completed_sourmash)
for filename in trim_files:
    #if filename.endswith("P.fq"):
    #if filename.endswith(".head"):
    if filename.endswith(".abundtrim"):
        filename_info = filename.split("_")
        genus_species_sample = filename_info[0] + "_" + filename_info[1] + "_" + filename_info[2] + "_" + filename_info[3][:1] + "_" + filename_info[4]
        trim_filename = trim_dir + filename
        sourmash_filename = trim_dir  + filename
        #head_reads(sourmash_filename,sourmash_dir,filename)
        #sourmash_files = [s for s in complete_sourmash_files if s.startswith(genus_species_sample)]
        #run_trim_low_abund(trim_filename,trim_low_abund_dir)
        #if len(sourmash_files) <1:
            #get_sourmash_command(sourmash_filename,sourmash_dir,genus_species_sample)
        get_sourmash_streaming(sourmash_filename,sourmash_dir,genus_species_sample)
