# Lisa Johnson
# commonly-used functions
# to use, "import clusterfunc"

import os
import subprocess
from subprocess import Popen, PIPE


def check_dir(dirname):
    if os.path.isdir(dirname) == False:
        os.mkdir(dirname)
        print("Directory created:", dirname)

def get_sbatch_filename(basedir,process_name,filename):
    sbatch_dir = basedir + "sbatch_files/"
    check_dir(sbatch_dir)
    sbatch_filename = sbatch_dir + process_name + "_" + filename + ".qsub"
    return sbatch_dir, sbatch_filename

def get_qsub_filename(basedir, process_name, filename):
    qsub_dir = basedir + "qsub_files/"
    check_dir(qsub_dir)
    qsub_filename = qsub_dir + process_name + "_" + filename + ".qsub"
    return qsub_dir, qsub_filename


def get_module_load_list(module_name_list):
    module_list = []
    for module in module_name_list:
        module_load = "module load " + module
        module_list.append(module_load)
    return module_list


def qsub_file(basedir, process_name, module_name_list, filename, process_string):
    working_dir = os.getcwd()
    qsub_dir, qsub_filename = get_qsub_filename(
        basedir, process_name, filename)
    os.chdir(qsub_dir)
    module_list = get_module_load_list(module_name_list)
# directive to use Laconia
# #PBS -l feature=intel16
#export MKL_NUM_THREADS=8
#export OMP_NUM_THREADS=8
# for ged lab nodes only
#PBS -A ged
    f = """#!/bin/bash
#PBS -l walltime=24:00:00,nodes=1:ppn=24
#PBS -l mem=500gb
#PBS -A ged
#PBS -j oe
cd ${{PBS_O_WORKDIR}}
""".format()
    with open(qsub_filename, "w") as qsub:
        qsub.write(f)
        for module_string in module_list:
            qsub.write(module_string + "\n")
        for string in process_string:
            qsub.write(string + "\n")
        qsub.write("qstat -f ${PBS_JOBID}\n")
        qsub.write(
            "cat ${PBS_NODEFILE} # Output Contents of the PBS NODEFILE\n")
        qsub.write(
            "env | grep PBS # Print out values of the current jobs PBS environment variables\n")
    qsub_string = 'qsub -V ' + qsub_filename
    print(qsub_string)
    #s = subprocess.Popen(qsub_string, shell=True)
    #s.wait()
    os.chdir(working_dir)                          

def sbatch_file(basedir,process_name,module_name_list,filename,process_string):
    working_dir=os.getcwd()
    sbatch_dir,sbatch_filename=get_sbatch_filename(basedir,process_name,filename)
    os.chdir(sbatch_dir)
    module_load=get_module_load_list(module_name_list)
    with open(sbatch_filename,"w") as sbatch_file:
        sbatch_file.write("#!/bin/bash -l"+"\n")
        sbatch_file.write("#SBATCH -D "+sbatch_dir+"\n")
        sbatch_file.write("#SBATCH -J "+process_name+"\n")
        sbatch_file.write("#SBATCH -o "+sbatch_dir+process_name+"-%j.o"+"\n")
        sbatch_file.write("#SBATCH -e "+sbatch_dir+process_name+"-%j.o"+"\n")
        sbatch_file.write("#SBATCH -t 240:00:00"+"\n")
        sbatch_file.write("#SBATCH -N 1"+"\n")
        #sbatch_file.write("#SBATCH -p RM"+"\n")
        sbatch_file.write("#SBATCH -p LM"+"\n")
        sbatch_file.write("#SBATCH -n 1"+"\n")
        sbatch_file.write("#SBATCH -c 48"+"\n")
        sbatch_file.write("#SBATCH --mem=2000GB"+"\n")
        #sbatch_file.write("#SBATCH --mem=4300MB"+"\n")
        for module_string in module_load:
            sbatch_file.write(module_string+"\n")
            #print(module_string)
        for string in process_string:
            sbatch_file.write(string+"\n")
            print(string)
    sbatch_string="sbatch --get-user-env "+sbatch_filename
    print(sbatch_string)
    s=subprocess.Popen(sbatch_string,shell=True)
    s.wait()
    os.chdir(working_dir)
