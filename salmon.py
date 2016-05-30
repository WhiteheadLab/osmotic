import os
import os.path
import subprocess
from subprocess import Popen, PIPE
# custom Lisa module
import clusterfunc

def get_data(thefile):
    count=0
    url_data={}
    with open(thefile,"rU") as inputfile:
        headerline=next(inputfile).split(',')
        #print headerline        
        position_name=headerline.index("ScientificName")
        position_reads=headerline.index("Run")
        position_ftp=headerline.index("download_path")
        for line in inputfile:
            line_data=line.split(',')
            name="_".join(line_data[position_name].split())
            read_type=line_data[position_reads]
            ftp=line_data[position_ftp]
            name_read_tuple=(name,read_type)
            print name_read_tuple
            #check to see if Scientific Name and run exist
            if name_read_tuple in url_data.keys():
                #check to see if ftp exists
                if ftp in url_data[name_read_tuple]:
                    print "url already exists:", ftp
                else:
                    url_data[name_read_tuple].append(ftp)
            else:
                url_data[name_read_tuple] = [ftp]
        return url_data

def salmon_index(newdir,genus_species,trinity_fasta):
	index=genus_species+"_index"
	salmon_index_string="salmon index --index "+newdir+index+" --transcripts "+trinity_fasta+" --type quasi"
	return salmon_index_string,index

def quant_salmon(newdir,dirname,genus_species,trinity_fasta,species):
	salmon_index_string,index=salmon_index(newdir,genus_species,trinity_fasta)
	print salmon_index_string
	salmon_string="""
for i in {}{}*.trim_1P.fq
do
	BASE=$(basename $i .trim_1P.fq)
	salmon quant -i {}{} --libType IU -1 {}$BASE.trim_1P.fq -2 {}$BASE.trim_2P.fq -o {}$BASE.quant;
done
""".format(dirname,species,newdir,index,dirname,dirname,newdir)
	print salmon_string
	salmonstring=[salmon_index_string,salmon_string]
        process_name="salmon"
        module_name_list=""
        clusterfunc.sbatch_file(newdir,process_name,module_name_list,genus_species,salmonstring)

	
def execute(assemblydirs,salmondir,assemblydir,basedir,trimdir):
	for genus_species_names in assemblydirs:
		genus_species = genus_species_names.split(".")[0]
		species = genus_species+genus_species_names.split(".")[1]
		print genus_species
		print species
		dirname=trimdir+genus_species+"/"
		newdir=salmondir+genus_species_names+"/"
		clusterfunc.check_dir(newdir)
		
		trinity_fasta=assemblydir+genus_species_names+"/"+genus_species_names+".Trinity.fixed.fa"
		quant_salmon(newdir,dirname,genus_species_names,trinity_fasta,species)


basedir="/home/ljcohen/osmotic/"
trimdir="/home/ljcohen/osmotic_trim_"
assemblydir="/home/ljcohen/msu_assemblies_finished/"
salmondir="/home/ljcohen/osmotic_salmon/"
clusterfunc.check_dir(trimdir)
clusterfunc.check_dir(assemblydir)
clusterfunc.check_dir(salmondir)
#assemblydirs=os.listdir(assemblydir)
assemblydirs=["F_heteroclitus.MDPL","F_heteroclitus.MDPP"]
execute(assemblydirs,salmondir,assemblydir,basedir,trimdir)
