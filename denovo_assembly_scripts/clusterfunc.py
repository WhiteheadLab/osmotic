# Lisa Cohen
# commonly-used functions
# to use, "import clusterfunc"

import os
import subprocess
from subprocess import Popen, PIPE

def check_dir(dirname):
    if os.path.isdir(dirname)==False:
        os.mkdir(dirname)
        print "Directory created:",dirname

def get_sbatch_filename(basedir,process_name,filename):
    sbatch_dir=basedir+"sbatch_files/"
    #sbatch_dir="/home/ljcohen/bin/"
    check_dir(sbatch_dir)
    sbatch_filename=sbatch_dir+process_name+"_"+filename+".sbatch"
    return sbatch_dir,sbatch_filename

def get_module_load_list(module_name_list):
    module_list=[]
    for module in module_name_list:
        module_load="module load "+module
        module_list.append(module_load)
    return module_list

def sbatch_file(basedir,process_name,module_name_list,filename,process_string):
    working_dir=os.getcwd()
    sbatch_dir,sbatch_filename=get_sbatch_filename(basedir,process_name,filename)
    os.chdir(sbatch_dir)
    module_load=get_module_load_list(module_name_list)
    with open(sbatch_filename,"w") as sbatch_file:
         sbatch_file.write("#!/bin/bash -l"+"\n")
         sbatch_file.write("#SBATCH -D "+sbatch_dir+"\n")
	 sbatch_file.write("#SBATCH -J "+process_name+"\n")
	 sbatch_file.write("#SBATCH -p high"+"\n")
	 sbatch_file.write("#SBATCH -t 04:00:00"+"\n")
	 #sbatch_file.write("#SBATHC -p med\n")
	 sbatch_file.write("#SBATCH --mail-user=ljcohen@ucdavis.edu"+"\n")
         sbatch_file.write("#SBATCH --cpus-per-task=16"+"\n")
	 sbatch_file.write("#SBATCH --mem=48GB"+"\n")
         for module_string in module_load:
             sbatch_file.write(module_string+"\n")
             #print module_string
         for string in process_string:
             sbatch_file.write(string+"\n")
             print string
    sbatch_string="sbatch "+sbatch_filename
    print sbatch_string
    s=subprocess.Popen(sbatch_string,shell=True)
    s.wait()
    os.chdir(working_dir)
