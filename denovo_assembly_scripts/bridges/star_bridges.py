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
                        else:
                                files_dictionary[sample_name]=[basedir+basefilename]
                        files_dictionary[sample_name] = sorted(files_dictionary[sample_name])
        return files_dictionary

def run_star(interleavedir,diginormdir,species):
        star_string="""
STAR --genomeDir /pylon5/bi5fpmp/ljcohen/Fhet_reference_genome/ensembl/Fhet_ensembl_star/ \
--runThreadN 8 \
--readFilesIn {}{} {}{} \
--outFileNamePrefix $species \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard 
""".format(interleavedir,species,diginormdir,species,diginormdir,species)
        print(stream_string)
        stream_command = [stream_string]
        module_load_list = []
        processname = "diginorm_stream"
        clusterfunc.sbatch_file(diginormdir,processname,module_load_list,species,stream_command)

def get_extract_paired(genus_species_dir,sample,abund_filt_filename):
	j = """
extract-paired-reads.py -p {}{}.pe.orphans -s {}{}.se.orphans {}{}.orphans.abundfilt
""".format(genus_species_dir,sample,genus_species_dir,sample,genus_species_dir,sample)
#        j="""
#extract-paired-reads.py -p {}{}.pe.keep.abundfilt.fq -s {}{}.se.keep.abundfilt.fq {}{}.orphans.abundfilt
#""".format(genus_species_dir,sample,genus_species_dir,sample,genus_species_dir,sample)
	return j

def run_extract_paired(genus_species_dir,sample):
	abund_filt_filename = genus_species_dir + sample + ".abundfilt"
        extract_paired=get_extract_paired(genus_species_dir,sample,abund_filt_filename)
        extract_command=[extract_paired]
        process_name="extract"
        module_name_list=""
        filename=sample
        clusterfunc.sbatch_file(genus_species_dir,process_name,module_name_list,filename,extract_command)

def get_genus_species(abund_filt_files):
	genus_species_files={}
	for filename in abund_filt_files:
		sample_info=filename.split("_")
		genus_species=sample_info[0]+"_"+sample_info[1]
		if genus_species in genus_species_files.keys():
			genus_species_files[genus_species].append(filename)
                        genus_species_files[genus_species] = sorted(genus_species_files[genus_species])
		else:
			genus_species_files[genus_species]=[filename]
	return genus_species_files

def split_paired(assemblyfilesdir,diginormdir,species):
	j = """ 
split-paired-reads.py -d {} {}{}.paired.gz
""".format(assemblyfilesdir,diginormdir,species)
	return j

def run_split_paired(assemblyfilesdir,diginormdir,species):
	split_command=split_paired(assemblyfilesdir,diginormdir,species)
        print(split_command)
        split_paired_command=[split_command]
        process_name="split"
        module_name_list=""
        filename=species
        clusterfunc.sbatch_file(diginormdir,process_name,module_name_list,filename,split_paired_command)

def combine_left(assemblyfilesdir,assemblydir,split_file,species):
	combine="""
cat {}{} > {}{}.left.fq
gunzip -c {}{}.single.gz >> {}{}.left.fq
gunzip -c {}{}.orphans.fq.gz >> {}{}.left.fq
""".format(assemblyfilesdir,split_file,assemblydir,species,assemblyfilesdir,species,assemblydir,species,assemblyfilesdir,species,assemblydir,species)
        combine_command=[combine]
	process_name="combine_left"
        module_name_list=""
	clusterfunc.sbatch_file(assemblyfilesdir,process_name,module_name_list,species,combine_command)

def combine_right(assemblyfilesdir,assemblydir,split_file,species):
        combine="""
cat {}{} > {}{}.right.fq
""".format(assemblyfilesdir,split_file,assemblydir,species)
        combine_command=[combine]
        process_name="combine_right"
        module_name_list=""
        clusterfunc.sbatch_file(assemblyfilesdir,process_name,module_name_list,species,combine_command)
def execute(listoffiles,trimdir,interleavedir,diginormdir,assemblyfilesdir,assemblydir):
        #files_dictionary=get_files(listoffiles,trimdir)
        #for sample in files_dictionary.keys():
        #        #print(sample)
        #        fileslist=sorted(files_dictionary[sample])
                #print(fileslist)
                #interleave_reads(trimdir,interleavedir,fileslist,sample)        
        interleavefileslist = os.listdir(interleavedir)
        #print(interleavefileslist)
        genus_species_files = get_genus_species(interleavefileslist)
        #print(genus_species_files)
        for species in genus_species_files:
            print(species)
            #orphansfile = combine_orphans(trimdir,sample)
            #streaming_diginorm(interleavedir,diginormdir,species)
        diginorm_files = os.listdir(diginormdir)
        for diginorm_file in diginorm_files:
            if diginorm_file.endswith(".paired.gz"):
                species = diginorm_file.split(".")[0]
                #run_split_paired(assemblyfilesdir,diginormdir,species) 
        split_files = os.listdir(assemblyfilesdir)
        for split_file in split_files:
            if split_file.endswith(".gz.1"):
                species = split_file.split(".")[0]
                combine_left(assemblyfilesdir,assemblydir,split_file,species)
            if split_file.endswith(".gz.2"):
                species = split_file.split(".")[0]
                combine_right(assemblyfilesdir,assemblydir,split_file,species) 
 
trimdir="/pylon5/bi5fpmp/ljcohen/kfish_trimmed/"
interleavedir="/pylon5/bi5fpmp/ljcohen/kfish_interleaved/"
diginormdir="/pylon5/bi5fpmp/ljcohen/kfish_diginorm/"
assemblyfilesdir="/pylon5/bi5fpmp/ljcohen/kfish_assemblies/"
assemblydir="/pylon5/bi5fpmp/ljcohen/kfish_txomes/"
clusterfunc.check_dir(interleavedir)
clusterfunc.check_dir(diginormdir)
clusterfunc.check_dir(assemblyfilesdir)
clusterfunc.check_dir(assemblydir)
listoffiles=os.listdir(trimdir)
execute(listoffiles,trimdir,interleavedir,diginormdir,assemblyfilesdir,assemblydir)
#group_assembly_files(diginormdir,assemblydir)
#combine_orphans(assemblydir)
#split_reads(assemblydir)
