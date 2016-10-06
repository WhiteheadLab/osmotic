import os
import os.path
from os.path import basename
from urllib import urlopen
from urlparse import urlparse
import subprocess
from subprocess import Popen, PIPE
import urllib
import shutil
import glob
import yaml
# custom Lisa module
import clusterfunc


def get_assemblies(assemblydir,trimdir):
	genus_species_dirs = os.listdir(assemblydir)
	trimfiles = os.listdir(trimdir)
	#genus_species_dirs=["F_heteroclitus.MDPP","F_heteroclitus.MDPL"]
	for genus_species in genus_species_dirs:
		print genus_species
		genus_species_trimfiles = get_reads(genus_species,trimfiles,trimdir)
		genus_species_dir=assemblydir+genus_species+"/"
		assemblyfile=genus_species_dir+genus_species+".Trinity.fixed.fa"
		if os.path.isfile(assemblyfile):
			print assemblyfile
		else:
			print "Not found:",assemblyfile
		print genus_species_trimfiles

# use diginorm reads to make it go faster
# too many reads for trimmed reads only
def get_reads(genus_species,trimfiles,trimdir):
	genus_species_trimfiles = []
	for filename in trimfiles:
		if filename.startswith(genus_species):
			genus_species_trimfiles.append(filename)
	return genus_species_trimfiles

def get_quant_file(salmondir, sra):
    quant_file = salmondir + sra + ".quant"
    return quant_file


def get_config_file(quant_file, rapclustdir, sra):
    config_file = rapclustdir + sra + "_config.yaml"
    out_dir = rapclustdir + sra + "_rapclust_out"
    yaml_dict = {"conditions": [sra], "samples": {
        sra: [quant_file]}, "outdir": out_dir}
    yaml_dump = yaml.dump(yaml_dict)
    with open(config_file, "w") as config_file:
        config_file.write(yaml_dump)
    return config_file


def run_rap_clust(salmondir, rapclustdir, sra):
    quant_file = get_quant_file(salmondir, sra)
    config_file = get_config_file(quant_file, rapclustdir, sra)
    config_filename = rapclustdir + sra + "_config.yaml"
    rapclust_string = "RapClust --config " + str(config_filename)
    print rapclust_string
    commands = [rapclust_string]
    process_name = "rapclust"
    module_name_list = ""
    filename = sra
    clusterfunc.qsub_file(rapclustdir, process_name,
                          module_name_list, filename, commands)


def execute(url_data):
    for item in url_data.keys():
        organism = item[0]
        org_seq_dir = basedir + organism + "/"
        url_list = url_data[item]
        for url in url_list:
            sra = basename(urlparse(url).path)
            newdir = org_seq_dir + sra + "/"
            salmondir = newdir + "salmon/"
            rapclustdir = newdir + "rapclust/"
            clusterfunc.check_dir(rapclustdir)
            clusterfunc.check_dir(salmondir)
            run_rap_clust(salmondir, rapclustdir, sra)

assemblydir = "/home/ljcohen/msu_assemblies_finished/"
trimdir = "/home/ljcohen/osmotic_trim/"
get_assemblies(assemblydir,trimdir)
