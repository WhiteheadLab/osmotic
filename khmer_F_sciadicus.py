import os
import os.path
from os.path import basename
import subprocess
from subprocess import Popen, PIPE
import clusterfunc

def parse_filename(filename):
        listoffilestomerge=[]
        fields=filename.split("_")
        genus=fields[0]
        species=fields[1]
        population=fields[2]
        treatment=fields[3]
        sample=fields[4]
        read=fields[5]
        extension=fields[6]
        sample=(genus,species,population,treatment,sample)
        return sample,genus_species

def get_files(listoffiles,basedir):
        files_dictionary={}
        for basefilename in listoffiles:
                if basefilename.endswith(".fq"):
                        filename=basedir+basefilename
                        fields=basefilename.split("_")
                        sample_name_info=fields[:-1]
                        sample_name="_".join(sample_name_info)
                        if sample_name in files_dictionary.keys():
                                files_dictionary[sample_name].append(basedir+basefilename)
                        else:
                                files_dictionary[sample_name]=[basedir+basefilename]
        return files_dictionary

def get_orphans(files_list,trimdir,sample):
	orphans_list=[]
	for filename in files_list:
		#print filename
        	if filename.endswith("1U.fq"):
			orphans_list.append(filename)
		elif filename.endswith("2U.fq"):
			orphans_list.append(filename)
        orphans_string=[make_orphans(trimdir,orphans_list,sample)]
        process_name="orphans"
        module_name_list=""
        filename=sample
        clusterfunc.sbatch_file(trimdir,process_name,module_name_list,filename,orphans_string)

def make_orphans(trimdir,orphans_list,sample):
	sorted_orphanslist=sorted(orphans_list)
    	R1=sorted_orphanslist[0]
	R2=sorted_orphanslist[1]
	orphan_string="""
gzip -9c {} {} > {}{}.orphans.fq.gz
""".format(R1,R2,trimdir,sample)
	return orphan_string

# all orphans are separate
# combine together with
# gzip -9c *.orphans > all.orphans.fq.gz

def get_pairs(files_list,trimdir,sample):
	paired_list=[]
	for filename in files_list:
		if filename.endswith("1P.fq"):
			paired_list.append(filename)
		elif filename.endswith("2P.fq"):
			paired_list.append(filename)
	return paired_list

def get_interleave_string(paired_list,interleavefile):
        sorted_pairs=sorted(paired_list)
        R1=sorted_pairs[0]
        R2=sorted_pairs[1]
        interleave_string="""
interleave-reads.py {} {} | gzip -9c > {}           
""".format(R1,R2,interleavefile)
        return interleave_string

def interleave_reads(trimdir,interleavedir,files_list,sample):
	interleavefile=interleavedir+sample+".interleaved.fq.gz"
	paired_list=get_pairs(files_list,trimdir,sample)
	interleave_string=[get_interleave_string(paired_list,interleavefile)]
	process_name="interleave"
	module_name_list=""
	filename=sample
	#clusterfunc.sbatch_file(interleavedir,process_name,module_name_list,filename,interleave_string)
	return interleavefile

def get_diginorm_string(diginorm_keep_file,diginormdir,orphansfile,interleavedir,genus_species):
# this will annoyingly put *.keep files in sbatch_files
# fix this!
	j="""
normalize-by-median.py -p -k 20 -C 20 -M 4e9 \\
--savegraph {}{}.norm.C20k20.ct \\
-u {} \\
{}*interleaved.fq.gz
""".format(diginormdir,genus_species,orphansfile,interleavedir)
	return j
	
def run_diginorm(diginormdir,orphansfile,interleavefile,sample):
		diginorm_keep_file=diginormdir+sample+".keep"
		diginorm_string=get_diginorm_string(diginorm_keep_file,diginormdir,orphansfile,interleavefile,sample)
		diginorm_command=[diginorm_string]
		process_name="diginorm"
		module_name_list=""
		filename=sample
		#clusterfunc.sbatch_file(diginormdir,process_name,module_name_list,filename,diginorm_command)
		graph_count_filename=diginormdir+sample+".norm.C20k20.ct"
		return graph_count_filename

def get_filter_abund(abundfilt_filename,diginormdir,genus_species):
	#if glob.glob(diginormdir+"*keep.abundfilt*"):
	#	print "filter-abund.py already run"
	#else:
	j="""
filter-abund.py -V -Z 18 -o {}.abundfilt \\
{}{}.norm.C20k20.ct \\
{}*.keep
""".format(abundfilt_filename,diginormdir,genus_species,diginormdir)
	return j

def run_filt_abund(diginormdir,graph_count_filename,genus_species):
	#if diginormfile.endswith("orphans.fq.gz.keep"):
	#	abundfilt_filename=diginormdir+diginorm_sample+".orphans"
	#else:
	abundfilt_filename=diginormdir+genus_species
	abund_filt=get_filter_abund(abundfilt_filename,diginormdir,genus_species)
	abund_filt_command=[abund_filt]
        process_name="abundfilt"
        module_name_list=""
        filename=genus_species
        clusterfunc.sbatch_file(diginormdir,process_name,module_name_list,filename,abund_filt_command)
	return abundfilt_filename
	
def get_rename(genus_species_dir,sample,abund_filt_filename):
        j="""
extract-paired-reads.py -p {}{}.pe.keep.abundfilt.fq -s {}{}.se.keep.abundfilt.fq {}
""".format(genus_species_dir,sample,genus_species_dir,sample,abund_filt_filename)
	return j

def get_diginorm_files(diginormdir):
        listoffiles=os.listdir(diginormdir)
        diginorm_files={}
        for filename in listoffiles:
                if filename.endswith(".keep"):
			file_info=filename.split("_")
			extension=file_info[-1].split(".")
			sample="_".join(file_info[:-1])
			sample=sample+"_"+extension[0]+".trim"
			if sample in diginorm_files.keys():
				diginorm_files[sample].append(filename)
			else:
				diginorm_files[sample]=[filename]
        return diginorm_files

def execute(listoffiles,trimdir,interleavedir,diginormdir,assemblydir):
        files_dictionary=get_files(listoffiles,trimdir)
	print files_dictionary
	for sample in files_dictionary.keys():
                fileslist=sorted(files_dictionary[sample])
		#get_orphans(fileslist,trimdir,sample)
		interleavefile=interleave_reads(trimdir,interleavedir,fileslist,sample)
# run diginorm for each genus_species
# you were running it per sample before
# does this make a difference?
	interleave_files=os.listdir(interleavedir)
	genus_species_files=get_genus_species(interleave_files)
	print genus_species_files
	for genus_species in genus_species_files: 
		print genus_species
		orphansfile=trimdir+genus_species+"_all.orphans.fq.gz"
		print os.path.isfile(orphansfile)
		print orphansfile
		print "If false, run this:"
		print "cd "+trimdir
		print "zcat *.orphans.fq.gz | gzip -9c > F_sciadicus_all.orphans.fq.gz"
		graph_count_filename=run_diginorm(diginormdir,orphansfile,interleavedir,genus_species)
	diginorm_files=get_diginorm_files(diginormdir)
	#print diginorm_files
	#for diginormfile in diginorm_files[sample]:
		#run_filt_abund(diginormdir,graph_count_filename,genus_species)
		
trimdir="/home/ljcohen/osmotic_trim_F_sciadicus/"
interleavedir=trimdir+"interleave/"
diginormdir=trimdir+"diginorm/"
assemblydir=trimdir+"assemblies/"
clusterfunc.check_dir(interleavedir)
clusterfunc.check_dir(diginormdir)
clusterfunc.check_dir(assemblydir)
listoffiles=os.listdir(trimdir)
#execute(listoffiles,trimdir,interleavedir,diginormdir,assemblydir)
#group_assembly_files(diginormdir,assemblydir)
#combine_orphans(assemblydir)
split_reads(assemblydir)
