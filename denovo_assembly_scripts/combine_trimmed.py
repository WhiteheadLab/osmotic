import os
import os.path
import subprocess
from subprocess import Popen,PIPE


def get_trimmed_files_left(basedir,genus_species,listoftrimmedfiles):
    left_files = []
    for filename in listoftrimmedfiles:
        if filename.startswith(genus_species):
            if filename.endswith("1P.fq"):
                left_files.append(basedir + filename)
    left_files.sort()
    return left_files


def get_trimmed_files_right(basedir,genus_species,listoftrimmedfiles):
    right_files = []
    for filename in listoftrimmedfiles:
        if filename.startswith(genus_species):
            if filename.endswith("2P.fq"):
                right_files.append(basedir + filename)
    right_files.sort()
    return right_files

def execute(assemblies_list,assemblies,listoftrimmedfiles,basedir,combined_dir):
    for assembly in assemblies_list:
        print(assembly)
        genus_species = assembly.split(".")[0]
        print(genus_species)
        trimmed_files_left = get_trimmed_files_left(basedir,genus_species,listoftrimmedfiles)
        trimmed_files_right = get_trimmed_files_right(basedir,genus_species,listoftrimmedfiles)
        out_file_left = combined_dir + genus_species + "_left.fq"
        out_file_right = combined_dir + genus_species + "_right.fq"
        files_left = " ".join(trimmed_files_left)
        files_right = " ".join(trimmed_files_right)
        command_left = "cat "+files_left + " > "+out_file_left
        command_right = "cat " +files_right+" > "+out_file_right
        print("Merging left...")
        print(command_left)
        s=subprocess.Popen(command_left,shell=True)
        s.wait()
        print("Merging right...")
        print(command_right)
        t=subprocess.Popen(command_right,shell=True)
        t.wait()

assemblies = "/mnt/home/ljcohen/osmotic_killifish_assemblies/assemblies/"
assemblies_list = os.listdir(assemblies)
basedir = "/mnt/scratch/ljcohen/osmotic_killifish/osmotic_trim/"
listoftrimmedfiles = os.listdir(basedir)
combined_dir = "/mnt/scratch/ljcohen/osmotic_killifish/combined_trimmed/"
execute(assemblies_list,assemblies,listoftrimmedfiles,basedir,combined_dir)
