import os
import os.path
from os.path import basename
import subprocess
from subprocess import Popen, PIPE
import clusterfunc

def get_files(listoffiles,basedir):
    files_dictionary={}
    for basefilename in listoffiles:
        if basefilename.endswith("P.fq"):
            filename=basedir+basefilename
            fields=basefilename.split("_")
            sample_name_info=fields[:-1]
            sample_name="_".join(sample_name_info)
            sample_name_split = sample_name.split(".")
            sample_name = sample_name_split[0]
            print(sample_name)
            if sample_name in files_dictionary.keys():
                files_dictionary[sample_name].append(basedir+basefilename)
                files_dictionary[sample_name] = sorted(files_dictionary[sample_name])
            else:
                files_dictionary[sample_name]=[basedir+basefilename]
    return files_dictionary

def run_star(sample,read1,read2,stardir):
        star_string="""
source /home/ljcohen/.bashrc
module load samtools
module load bwa
cd /pylon5/bi5fpmp/ljcohen/kfish_bwa
bwa mem \\
/pylon5/bi5fpmp/ljcohen/Fhet_reference_genome/ncbi/Fhet_ncbi_bwa/GCF_000826765.1_Fundulus_heteroclitus-3.0.2_genomic.fna \\
{} \\
{} \\
> {}{}.sam
samtools sort -@ 8 -o {}{}.bam {}{}.sam
rm {}{}.sam
samtools flagstat {}{}.bam > {}{}.bam.stats
""".format(read1,read2,stardir,sample,stardir,sample,stardir,sample,stardir,sample,stardir,sample,stardir,sample)
        print(star_string)
        star_command = [star_string]
        module_load_list = []
        processname = "bwa_"+sample
        clusterfunc.sbatch_file(stardir,processname,module_load_list,sample,star_command)

def execute(listoffiles,trimdir,stardir):
        files_dictionary=get_files(listoffiles,trimdir)
        print(files_dictionary)
        for sample in files_dictionary.keys():
            read1=files_dictionary[sample][0]
            read2=files_dictionary[sample][1]
            run_star(sample,read1,read2,stardir)
 
trimdir="/pylon5/bi5fpmp/ljcohen/kfish_trimmed/"
stardir="/pylon5/bi5fpmp/ljcohen/kfish_bwa/"
listoffiles=os.listdir(trimdir)
execute(listoffiles,trimdir,stardir)
