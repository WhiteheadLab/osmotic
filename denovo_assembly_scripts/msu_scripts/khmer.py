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
		print filename
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
/home/ljcohen/bin/khmer/scripts/interleave-reads.py {} {} | gzip > {}           
""".format(R1,R2,interleavefile)
        return interleave_string

def interleave_reads(trimdir,interleavedir,files_list,sample):
	interleavefile=interleavedir+sample+".interleaved.fq.gz"
	paired_list=get_pairs(files_list,trimdir,sample)
	interleave_string=[get_interleave_string(paired_list,interleavefile)]
	process_name="interleave"
	module_name_list=""
	filename=sample
	#clusterfunc.sbatch_file(trimdir,process_name,module_name_list,filename,interleave_string)
	return interleavefile

def get_diginorm_string(diginorm_keep_file,diginormdir,orphansfile,interleavefile,sample):
	j="""
normalize-by-median.py -p -k 20 -C 20 -M 4e9 \\
--savegraph {}{}.norm.C20k20.ct \\
-u {} \\
{}
""".format(diginormdir,sample,orphansfile,interleavefile)
	return j
	
def run_diginorm(diginormdir,orphansfile,interleavefile,sample):
		diginorm_keep_file=diginormdir+sample+".keep"
		diginorm_string=get_diginorm_string(diginorm_keep_file,diginormdir,orphansfile,interleavefile,sample)
		diginorm_command=[diginorm_string]
		process_name="diginorm2"
		module_name_list=""
		filename=sample
		#clusterfunc.sbatch_file(diginormdir,process_name,module_name_list,filename,diginorm_command)
		graph_count_filename=diginormdir+sample+".norm.C2-k2-.ct"
		return graph_count_filename

def get_filter_abund(abundfilt_filename,diginormfile,diginormdir,sample):
	#if glob.glob(diginormdir+"*keep.abundfilt*"):
	#	print "filter-abund.py already run"
	#else:
	j="""
filter-abund.py -V -Z 18 -o {}.abundfilt \\
{}{}.norm.C20k20.ct \\
{}{}
""".format(abundfilt_filename,diginormdir,sample,diginormdir,diginormfile)
	return j

def run_filt_abund(diginormdir,graph_count_filename,diginormfile,diginorm_sample):
	if diginormfile.endswith("orphans.fq.gz.keep"):
		abundfilt_filename=diginormdir+diginorm_sample+".orphans"
	else:
		abundfilt_filename=diginormdir+diginorm_sample
	abund_filt=get_filter_abund(abundfilt_filename,diginormfile,diginormdir,diginorm_sample)
	abund_filt_command=[abund_filt]
        process_name="abundfilt"
        module_name_list=""
        filename=diginorm_sample
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
	diginorm_files=get_diginorm_files(diginormdir)
	for sample in files_dictionary.keys():
                fileslist=sorted(files_dictionary[sample])
		#get_orphans(fileslist,trimdir,sample)
		orphansfile=trimdir+sample+".orphans.fq.gz"
		interleavefile=interleave_reads(trimdir,interleavedir,fileslist,sample)
		graph_count_filename=run_diginorm(diginormdir,orphansfile,interleavefile,sample)
		diginorm_files=get_diginorm_files(diginormdir)
		for diginormfile in diginorm_files[sample]:
			run_filt_abund(diginormdir,graph_count_filename,diginormfile,sample)
		

def get_abund_filt_files(diginormdir):
	listoffiles=os.listdir(diginormdir)
	abund_filt_files=[]
	for filename in listoffiles:
		if filename.endswith(".abundfilt"):
			abund_filt_files.append(filename)
	return abund_filt_files

def get_genus_species(abund_filt_files):
	genus_species_files={}
	for filename in abund_filt_files:
		sample_info=filename.split("_")
		genus_species=sample_info[0]+"_"+sample_info[1]
		if genus_species in genus_species_files.keys():
			genus_species_files[genus_species].append(filename)
		else:
			genus_species_files[genus_species]=[filename]
	return genus_species_files

def consolidate(assemblydir,genus_species):
	j="""
echo -n "" > {}{}.orphans.keep.abundfilt.fq.gz

for file1 in {}*orphans.se.keep.abundfilt.fq
do
	gzip -9c ${{file1}} >> {}{}.orphans.keep.abundfilt.fq.gz
done

for file2 in {}*abundfilt.se.keep.abundfilt.fq
do
	gzip -9c ${{file2}} >> {}{}.orphans.keep.abundfilt.fq.gz
done
""".format(assemblydir,genus_species,assemblydir,assemblydir,genus_species,assemblydir,assemblydir,genus_species)
	return j

def group_assembly_files(diginormdir,assemblydir):		
	abund_filt_files=get_abund_filt_files(diginormdir)
	genus_species_files=get_genus_species(abund_filt_files)
	for genus_species in genus_species_files:
		genus_species_dir=assemblydir+genus_species+"/"
		clusterfunc.check_dir(genus_species_dir)
		for filename in genus_species_files[genus_species]:
			abund_filt_filename=diginormdir+filename
			sample_info=filename.split("_")
			extension=sample_info[-1].split(".")
			sample="_".join(sample_info[:-1])
			sample=sample+"_"+extension[0]+"_"+extension[2]
			rename_command=get_rename(genus_species_dir,sample,abund_filt_filename)
			process_name="split"
			rename_command=[rename_command]
        		module_name_list=""
        		filename=sample
        		#clusterfunc.sbatch_file(genus_species_dir,process_name,module_name_list,filename,rename_command)
		
def combine_orphans(assemblydir):
	assemblydirs=os.listdir(assemblydir)
	for genus_species in assemblydirs:
		genus_species_dir=assemblydir+genus_species+"/"
		listoffiles=os.listdir(genus_species_dir)
		consolidate_command=consolidate(genus_species_dir,genus_species)			
		consolidate_command=[consolidate_command]
		process_name="gzip"
		module_name_list=""
		filename=genus_species
		#clusterfunc.sbatch_file(genus_species_dir,process_name,module_name_list,filename,consolidate_command)
		for filename in listoffiles:
			if filename.endswith("pe.keep.abundfilt.fq"):
				gzip_command="""
gzip {}{}
""".format(genus_species_dir,filename)
				print gzip_command
				#s=subprocess.Popen(gzip_command,shell=True)
				#s.wait()

def split_reads(assemblydir):
	assemblydirs=os.listdir(assemblydir)
	for genus_species in assemblydirs:
		genus_species_dir=assemblydir+genus_species+"/"
		listoffiles=os.listdir(genus_species_dir)
		for filename in listoffiles:
			if filename.endswith("pe.keep.abundfilt.fq.gz"):
# next time you run this, specify output file,
# otherwise output will be put in sbatch_files directory
# moved manually 2/7/2016
				split_command="""
split-paired-reads.py {}{}
""".format(genus_species_dir,filename)
				process_name="split"
				module_name_list=""
				split_command=[split_command]
				#clusterfunc.sbatch_file(genus_species_dir,process_name,module_name_list,filename,split_command)
		combine="""
cat {}*.1 > {}{}.left.fq
cat {}*.2 > {}{}.right.fq
gunzip -c {}*orphans.keep.abundfilt.fq.gz >> {}{}.left.fq
""".format(genus_species_dir,genus_species_dir,genus_species,genus_species_dir,genus_species_dir,genus_species,genus_species_dir,genus_species_dir,genus_species)
		combine_command=[combine]
		process_name="combine"
                module_name_list=""
		clusterfunc.sbatch_file(genus_species_dir,process_name,module_name_list,genus_species,combine_command)

trimdir="/home/ljcohen/osmotic_trim/"
interleavedir="/home/ljcohen/osmotic_trim/interleave/"
diginormdir="/home/ljcohen/osmotic_trim/diginorm/"
assemblydir="/home/ljcohen/osmotic_trim/assemblies/"
clusterfunc.check_dir(interleavedir)
clusterfunc.check_dir(diginormdir)
clusterfunc.check_dir(assemblydir)
listoffiles=os.listdir(trimdir)
#execute(listoffiles,trimdir,interleavedir,diginormdir,assemblydir)
#group_assembly_files(diginormdir,assemblydir)
#combine_orphans(assemblydir)
split_reads(assemblydir)
