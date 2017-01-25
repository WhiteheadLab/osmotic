import os
import os.path
from os.path import basename
import subprocess
from subprocess import Popen, PIPE
import shutil
import glob
import yaml
from yaml.representer import Representer
# custom Lisa module
import clusterfunc_py3

def get_samples_dict(genus_species,salmondir):
    print(genus_species)
    salmonfiles = os.listdir(salmondir)
    genus_species_quant_dict = {}
    for quant_file_short in salmonfiles:
        if quant_file_short.endswith(".quant"):
                info = quant_file_short.split("_")
                #print(info)
                genus_species_2 = info[0] + "_" + info[1]
                if genus_species_2 == genus_species:
                    condition = info[2]
                    sample = info[3]
                    quant_file = salmondir + quant_file_short
                    if condition in genus_species_quant_dict:
                        #print(condition)
                        if quant_file in genus_species_quant_dict[condition]:
                            print("quant_file already exists:",quant_file)
                        else:
                            genus_species_quant_dict[condition].append(quant_file)
                    else:
                        #print(condition)
                        genus_species_quant_dict[condition]=[quant_file]
    #print(genus_species_quant_dict)
    return genus_species_quant_dict 

def get_config_file(quant_files, rapclustdir, genus_species):
    config_file = rapclustdir + genus_species + "_config.yaml"
    out_dir = rapclustdir + genus_species + "_rapclust_out"
    conditions = list(quant_files.keys())
    num_conditions = len(conditions)
    yaml_dict = {"conditions": conditions, "outdir": out_dir}
    yaml_dict["samples"] = quant_files
    print(yaml_dict)
    yaml_dump = yaml.dump(yaml_dict)
    with open(config_file, "w") as config_file:
        config_file.write(yaml_dump)
    return config_file

def run_rap_clust(salmondir, rapclustdir, genus_species):
    quant_files = get_samples_dict(genus_species,salmondir)
    config_file = get_config_file(quant_files, rapclustdir, genus_species)
    config_filename = rapclustdir + genus_species + "_config.yaml"
    rapclust_string = "RapClust --config " + str(config_filename)
    env_string = "source activate rapclust"
    commands = [env_string,rapclust_string]
    process_name = "rapclust"
    module_name_list = ""
    filename = genus_species
    #clusterfunc_py3.sbatch_file(rapclustdir, process_name,module_name_list, filename, commands)
    print(rapclust_string)
    s = subprocess.Popen(rapclust_string, shell=True)
    s.wait()

def execute(assemblies,trimdir,rapclustdir,salmondir):
    for assembly in assemblies:
        if assembly.endswith(".dammit.fasta"):
            genus_species = assembly.split(".")[0]
            #print(assembly)
            #print(genus_species)
            run_rap_clust(salmondir, rapclustdir, genus_species)

assemblydir = "/home/ljcohen/osmotic_assemblies_farm/"
rapclustdir = "/home/ljcohen/osmotic_rapclust/"
salmondir = "/home/ljcohen/osmotic_salmon/"
clusterfunc_py3.check_dir(rapclustdir)
trimdir = "/home/ljcohen/osmotic_trim/"
assemblies = os.listdir(assemblydir)
execute(assemblies,trimdir,rapclustdir,salmondir)
