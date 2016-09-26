import os
import os.path
from os.path import basename
import subprocess
from subprocess import Popen, PIPE
import clusterfunc

def fix_filename(trimdir,filename):
	if filename.startswith("F_heteroclitus"):
		print filename
		if filename.endswith(".fq") or filename.endswith(".fq.gz"):
        		if filename.endswith("P.fq") or filename.endswith ("U.fq"): 
				fields=filename.split(".")
        			sample_info=fields[0]
				trim=fields[1]
				extension1=fields[2]
				extension=trim+"."+extension1
        			fields2=sample_info.split("_")
				genus=fields2[0]
				species=fields2[1]
        			population=fields2[2]
        			species_population=species+"."+population
				treatment=fields2[3]
       	 			sample=fields2[4]
        			new_filename="_".join([genus,species_population,treatment,sample,extension])
			else:
				fields=filename.split("_")
				genus=fields[0]
				species=fields[1]
				population=fields[2]
				treatment=fields[3]
				extension=fields[4]
				species_population=species+"."+population
				new_filename="_".join([genus,species_population,treatment,extension])
        		fix_filename_command=""" 
mv {}{} {}{}
""".format(trimdir,filename,trimdir,new_filename)
			#s=subprocess.Popen(fix_filename_command,shell=True)
			#s.wait()
			print fix_filename_command

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
                        if basefilename.endswith("P.fq"):
				if sample_name in files_dictionary.keys():
                                	files_dictionary[sample_name].append(basedir+basefilename)
                        	else:
                                	files_dictionary[sample_name]=[basedir+basefilename]
        return files_dictionary


def get_files_population(listoffiles,basedir):
        files_dictionary={}
        for filename in listoffiles:
                if filename.startswith("F_heteroclitus"):
			if filename.endswith(".interleaved.fq.gz"):
				file_info=filename.split("_")
				population=file_info[0]+"_"+file_info[1]
                		full_filename=basedir+filename
                		if population in files_dictionary.keys():
                			files_dictionary[population].append(full_filename)
                		else:
                        		files_dictionary[population]=[full_filename]
        return files_dictionary


# orphans are all in ~/osmotic_trim
# ran this:
# zcat F_heteroclitus.MDPL*orphans* > interleave_fixed/F_heteroclitus.MDPL.trim.orphans.fq
# zcat F_heteroclitus.MDPP*orphans* > interleave_fixed/F_heteroclitus.MDPP.trim.orphans.fq

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
        #clusterfunc.sbatch_file(trimdir,process_name,module_name_list,filename,orphans_string)

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
	clusterfunc.sbatch_file(interleavedir,process_name,module_name_list,filename,interleave_string)
	return interleavefile

def get_diginorm_string(diginorm_keep_file,diginormdir,orphansfile,interleavefile,sample):
	j="""
normalize-by-median.py -p -k 20 -C 20 -M 4e9 \\
--savegraph {}{}.norm.C20k20.ct \\
-u {} \\
{}{}*.interleaved.fq.gz

""".format(diginormdir,sample,orphansfile,interleavedir,sample)
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

def get_filter_abund(diginormdir,sample):
	#if glob.glob(diginormdir+"*keep.abundfilt*"):
	#	print "filter-abund.py already run"
	#else:
	j="""
filter-abund.py -V -Z 18 -o {}{}.abundfilt \\
{}{}.norm.C20k20.ct \\
{}{}*.keep
""".format(diginormdir,sample,diginormdir,sample,diginormdir,sample)
	return j

def run_filt_abund(diginormdir,graph_count_filename,diginorm_sample):
#	if diginormfile.endswith("orphans.fq.gz.keep"):
#		abundfilt_filename=diginormdir+diginorm_sample+".orphans"
#	else:
#		abundfilt_filename=diginormdir+diginorm_sample
	abund_filt=get_filter_abund(diginormdir,diginorm_sample)
	abund_filt_command=[abund_filt]
	process_name="abundfilt"
        module_name_list=""
        filename=diginorm_sample
        clusterfunc.sbatch_file(diginormdir,process_name,module_name_list,filename,abund_filt_command)
	
def get_rename(diginormdir,sample):
        j="""
extract-paired-reads.py {}{}.abundfilt
""".format(diginormdir,sample)
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
        #files_dictionary=get_files(listoffiles,trimdir)
	#print files_dictionary
	#for sample in files_dictionary:
		#fileslist=sorted(files_dictionary[sample])
		#interleavefile=interleave_reads(trimdir,interleavedir,fileslist,sample)
        interleave_files=os.listdir(interleavedir)
	interleave_files_dictionary=get_files_population(interleave_files,interleavedir)	
	print interleave_files_dictionary
	#diginorm_files=get_diginorm_files(diginormdir)
	for population in interleave_files_dictionary.keys():
                #fileslist=sorted(files_dictionary[population])
		#get_orphans(fileslist,trimdir,sample)
		#orphansfile=interleavedir+population+".trim.orphans.fq.gz"
		#graph_count_filename=run_diginorm(diginormdir,orphansfile,interleavedir,population)
		#diginorm_files=get_diginorm_files(diginormdir)
		#for diginormfile in diginorm_files[sample]:
		#run_filt_abund(diginormdir,graph_count_filename,population)
		genus_species_dir=assemblydir+population+"/"
		rename_command=get_rename(diginormdir,population)
                process_name="split"
                rename_command=[rename_command]
                module_name_list=""
                filename=population
                
		#clusterfunc.sbatch_file(diginormdir,process_name,module_name_list,filename,rename_command)
		split_command1=split_reads(diginormdir,population)
		process_name="split"
        	module_name_list=""
        	split_command=[split_command1]
        	clusterfunc.sbatch_file(diginormdir,process_name,module_name_list,filename,split_command)


def split_reads(diginormdir,population):
	split_command="""
split-paired-reads.py -1 {}{}.1 -2 {}{}.2 {}{}.abundfilt.pe
""".format(diginormdir,population,diginormdir,population,diginormdir,population)
	return split_command	


def combine_split():
	combine="""
cat {}*.1 > {}{}.left.fq
cat {}*.2 > {}{}.right.fq
gunzip -c {}*orphans.keep.abundfilt.fq >> {}{}.left.fq
""".format(genus_species_dir,genus_species_dir,genus_species,genus_species_dir,genus_species_dir,genus_species,genus_species_dir,genus_species_dir,genus_species)
	combine_command=[combine]
	process_name="combine"
        module_name_list=""
	clusterfunc.sbatch_file(genus_species_dir,process_name,module_name_list,genus_species,combine_command)

trimdir="/home/ljcohen/osmotic_trim/"
interleavedir="/home/ljcohen/osmotic_trim/interleave_fixed/"
diginormdir="/home/ljcohen/osmotic_trim/diginorm_fixed/"
assemblydir="/home/ljcohen/osmotic_assemblies_completed/"
clusterfunc.check_dir(interleavedir)
clusterfunc.check_dir(diginormdir)
clusterfunc.check_dir(assemblydir)
listoffiles_all=os.listdir(trimdir)

#for filename in listoffiles:
#	new_filename=fix_filename(trimdir,filename)

listoffiles=[]
for filename in listoffiles_all:
	if filename.startswith("F_heteroclitus"):
		listoffiles.append(filename) 

#print listoffiles
	

execute(listoffiles,trimdir,interleavedir,diginormdir,assemblydir)
#group_assembly_files(diginormdir,assemblydir)
#combine_orphans(assemblydir)
#split_reads(assemblydir)
