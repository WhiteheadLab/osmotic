import os
import os.path
import subprocess
from subprocess import Popen, PIPE
import shutil
import glob
# custom Lisa module
import clusterfunc_py3


def fastqc_report(fastqcdir,fastq_file,sample_name):
    fastqc_string="fastqc -o "+fastqcdir+" "+fastq_file
    process_string=[fastqc_string]
    process_name="fastqc"
    module_load_list=["fastqc/0.10.1"]
    clusterfunc_py3.sbatch_file(fastqcdir,process_name,module_load_list,sample_name,process_string)

def execute(listoffiles,fastqcdir,trimdir):
    for filename in listoffiles:
        if filename.endswith("P.fq"):
            fastq_file = trimdir + filename
            fastqc_report(fastqcdir,fastq_file,filename)


trimdir="/home/ljcohen/osmotic_trim/"
fastqcdir="/home/ljcohen/osmotic_aftertrim_fastqc/"
listoffiles=os.listdir(trimdir)
execute(listoffiles,fastqcdir,trimdir)
