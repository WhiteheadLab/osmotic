import os
import os.path
from os.path import basename
import subprocess
from subprocess import Popen, PIPE
import clusterfunc

def get_files(listoffiles,basedir):
        files_dictionary={}
        for basefilename in listoffiles:
                if basefilename.endswith(".fq"):
                        filename=basedir+basefilename
                        fields=basefilename.split("_")
                        sample_name_info=fields[:-1]
                        sample_name="_".join(sample_name_info)
			#print sample_name
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
	print files_list
	paired_list=[]
	for filename in files_list:
		if filename.endswith("1P.fq"):
			paired_list.append(filename)
		elif filename.endswith("2P.fq"):
			paired_list.append(filename)
	return paired_list

def get_interleave_string(paired_list,interleavefile):
        print paired_list
	sorted_pairs=sorted(paired_list)
        R1=sorted_pairs[0]
        R2=sorted_pairs[1]
        interleave_string="""
interleave-reads.py {} {} | gzip > {}           
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

def get_diginorm_string(diginorm_keep_file,diginormdir,orphansfile,interleavedir,sample):
	j="""
normalize-by-median.py -p -k 20 -C 20 -M 4e9 \\
--savegraph {}{}.norm.C20k20.ct \\
-u {} \\
{}{}*.interleaved.fq.gz
""".format(diginormdir,sample,orphansfile,interleavedir,sample)
	return j
	
def run_diginorm(trimdir,diginormdir,orphansfile,interleaveadir,sample):
		diginorm_keep_file=diginormdir+sample+".keep"
		diginorm_string=get_diginorm_string(diginorm_keep_file,diginormdir,orphansfile,interleavedir,sample)
		diginorm_command=[diginorm_string]
		process_name="diginorm"
		module_name_list=""
		filename=sample
		#clusterfunc.sbatch_file(trimdir,process_name,module_name_list,filename,diginorm_command)
		graph_count_filename=diginormdir+sample+".norm.C20k20.ct"

def get_filter_abund(diginormdir,sample):
	#if glob.glob(diginormdir+"*keep.abundfilt*"):
	#	print "filter-abund.py already run"
	#else:
	j="""
filter-abund.py -V -Z 18 -o {}{}.orphans.abundfilt \\
{}{}.norm.C20k20.ct \\
{}*.orphans.fq.gz.keep
""".format(diginormdir,sample,diginormdir,sample,diginormdir)
	return j

def run_filt_abund(sample_diginormdir,sample):
	abund_filt=get_filter_abund(sample_diginormdir,sample)
	abund_filt_command=[abund_filt]
        process_name="abundfilt"
        module_name_list=""
        filename=sample
        clusterfunc.sbatch_file(trimdir,process_name,module_name_list,filename,abund_filt_command)
	
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

def get_samples(interleavefileslist,interleavedir):
	interleave_files = {}
	for filename in interleavefileslist:
		if filename.endswith(".fq.gz"):
			file_info = filename.split("_")
			if "NA" in file_info:
				position_NA = file_info.index("NA")
				file_info.pop(position_NA)
				sample_num = file_info[-1].split(".")[0]
				file_info.pop(-1)
				file_info.pop(-1)
				#file_info.append(sample_num)
				#print file_info
				sample = "_".join(file_info)
				#print sample
				
			else:
				sample_num = file_info[-2]
				file_info.pop(-1)
				file_info.pop(-1)
				file_info.pop(-1)
				#print file_info
				sample = "_".join(file_info)
				#print sample
			if sample in interleave_files:
				if filename in interleave_files[sample]:
					print "exists",filename
				else:
					interleave_files[sample].append(filename)
			else:
				interleave_files[sample] = [filename]
	return interleave_files			

def combine_orphans(trimdir,sample):
	orphansfile = trimdir + sample + ".orphans.fq.gz"
	j = """
touch {}
for file in {}{}*.orphans.fq
do
        gzip -9c ${{file}} >> {}
done 
""".format(orphansfile,trimdir,sample,orphansfile)
	#print j
	#s = subprocess.Popen(j, shell=True)
        #s.wait()	
	return orphansfile

def execute(listoffiles,trimdir,interleavedir,diginormdir,assemblydir):
        files_dictionary=get_files(listoffiles,trimdir)
	for sample in files_dictionary.keys():
                fileslist=sorted(files_dictionary[sample])
		matching = [s for s in fileslist if "orphans" in s]
	interleavefileslist = os.listdir(interleavedir)
	interleave_files = get_samples(interleavefileslist,interleavedir)
	print interleave_files
	for sample in interleave_files:
		orphansfile = combine_orphans(trimdir,sample)
		sample_diginormdir = diginormdir + sample + "/"
		clusterfunc.check_dir(sample_diginormdir)
		genus_species_assemblydir = assemblydir + sample + "/"
		clusterfunc.check_dir(genus_species_assemblydir)	
		#run_diginorm(trimdir,sample_diginormdir,orphansfile,interleavedir,sample)
		#run_filt_abund(sample_diginormdir,sample)
		#run_extract_paired(sample_diginormdir,sample)
		#combine_orphans_after_diginorm(sample_diginormdir,sample)
		#run_split_paired(sample_diginormdir,sample)
		combine(genus_species_assemblydir,sample_diginormdir,sample)

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

def consolidate(genus_species_diginormdir,sample):
	j="""
gzip -9c {}{}.orphans.abundfilt > {}{}.orphans.keep.abundfilt.fq.gz
gzip -9c {}{}.se.keep.abundfilt.fq >> {}{}.orphans.keep.abundfilt.fq.gz
""".format(genus_species_diginormdir,sample,genus_species_diginormdir,sample,genus_species_diginormdir,sample,genus_species_diginormdir,sample)
	return j


def split_paired(genus_species_diginormdir,sample):
	j = """ 
split-paired-reads.py -d {} {}{}.pe.keep.abundfilt.fq
""".format(genus_species_diginormdir,genus_species_diginormdir,sample)
	return j

def run_split_paired(genus_species_diginormdir,sample):
	split_command=split_paired(genus_species_diginormdir,sample)
        split_paired_command=[split_command]
        process_name="split"
        module_name_list=""
        filename=sample
        clusterfunc.sbatch_file(genus_species_diginormdir,process_name,module_name_list,filename,split_paired_command)

		
def combine_orphans_after_diginorm(genus_species_diginormdir,sample):
	consolidate_command=consolidate(genus_species_diginormdir,sample)			
	consolidate_command=[consolidate_command]
	process_name="gzip"
	module_name_list=""
	filename=sample
	clusterfunc.sbatch_file(genus_species_diginormdir,process_name,module_name_list,filename,consolidate_command)

def combine(genus_species_assemblydir,genus_species_diginormdir,sample):
		combine="""
cat {}{}*.1 > {}{}.left.fq
cat {}{}*.2 > {}{}.right.fq
gunzip -c {}{}.orphans.keep.abundfilt.fq.gz >> {}{}.left.fq
""".format(genus_species_diginormdir,sample,genus_species_assemblydir,sample,genus_species_diginormdir,sample,genus_species_assemblydir,sample,genus_species_diginormdir,sample,genus_species_assemblydir,sample)
		combine_command=[combine]
		process_name="combine"
                module_name_list=""
		clusterfunc.sbatch_file(genus_species_diginormdir,process_name,module_name_list,sample,combine_command)

trimdir="/home/ljcohen/osmotic_trim/"
interleavedir="/home/ljcohen/osmotic_trim/interleaved/"
diginormdir="/home/ljcohen/osmotic_trim/diginorm/"
assemblydir="/home/ljcohen/assemblies/"
clusterfunc.check_dir(interleavedir)
clusterfunc.check_dir(diginormdir)
clusterfunc.check_dir(assemblydir)
listoffiles=os.listdir(trimdir)
execute(listoffiles,trimdir,interleavedir,diginormdir,assemblydir)
#group_assembly_files(diginormdir,assemblydir)
#combine_orphans(assemblydir)
#split_reads(assemblydir)
