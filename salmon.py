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

def salmon_index(salmondir,genus_species,trinity_fasta,basedir):
	index=genus_species+"_index"
	salmondir_genus_species=salmondir+genus_species+"/"
	#os.chdir(salmondir_genus_species)
	salmon_index="salmon index --index "+index+" --transcripts "+trinity_fasta+" --type quasi"
	#print salmon_index	
	#s=subprocess.Popen(salmon_index,shell=True)
	#s.wait()
	#print "Indexed."
	#os.chdir(basedir)
	return index

def quant_salmon(newdir,index,genus_species,trimdir):
	#os.chdir(salmondir)
	salmon_string="""
for i in {}{}*.trim_1P.fq
do
	BASE=$(basename $i .trim_1P.fq)
	salmon quant -i {}{} --libType IU -1 {}$BASE.trim_1P.fq -2 {}$BASE.trim_2P.fq -o {}$BASE.quant;
done
""".format(trimdir,genus_species,newdir,index,trimdir,trimdir,newdir,genus_species)
	print salmon_string
        #s=subprocess.Popen(salmon_string,shell=True)
	#s.wait(i)
	salmonstring=[salmon_string]
        process_name="salmon"
        module_name_list=""
        clusterfunc.sbatch_file(newdir,process_name,module_name_list,genus_species,salmonstring)

	
def execute(assemblydirs,salmondir,assemblydir,basedir,trimdir):
	for genus_species in assemblydirs:
		newdir=salmondir+genus_species+"/"
		clusterfunc.check_dir(newdir)
		trinity_fasta=assemblydir+genus_species+"/"+genus_species+".Trinity.fixed.fa"
		index=salmon_index(salmondir,genus_species,trinity_fasta,basedir)
		quant_salmon(newdir,index,genus_species,trimdir)


basedir="/home/ljcohen/osmotic/"
trimdir="/home/ljcohen/osmotic_trim/"
assemblydir="/home/ljcohen/msu_assemblies_finished/"
salmondir="/home/ljcohen/osmotic_salmon/"
clusterfunc.check_dir(trimdir)
clusterfunc.check_dir(assemblydir)
clusterfunc.check_dir(salmondir)
assemblydirs=os.listdir(assemblydir)
execute(assemblydirs,salmondir,assemblydir,basedir,trimdir)
