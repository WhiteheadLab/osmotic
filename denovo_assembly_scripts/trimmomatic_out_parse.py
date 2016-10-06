import os

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
	return sample


def get_pairs(listoffiles,basedir):
	pairs_dictionary={}
	for basefilename in listoffiles:
		if basefilename.endswith(".fastq.gz"):
                	filename=basedir+basefilename
                	fields=basefilename.split("_")
                	sample_name_info=fields[:-1]
			sample_name="_".join(sample_name_info)
			if sample_name in pairs_dictionary.keys():
				pairs_dictionary[sample_name].append(basedir+basefilename)
			else:
				pairs_dictionary[sample_name]=[basedir+basefilename]
	return pairs_dictionary

def parse_trimmomatic(trim_out_file):
	important_nums = []
	if os.path.isfile(trim_out_file)==True:
		with open(trim_out_file) as outfile:
			for line_full in outfile:
				line=line_full.strip()
				if line.startswith("Input Read Pairs:") == True:
					reads_num=line.split()
					input_reads=reads_num[3]
					surviving_reads=reads_num[6]
					print input_reads
					important_nums.extend(input_reads, surviving_reads)	
				
	print important_nums
	return important_nums
				
def get_sample_dictionary(sample_dictionary,trim_out_file,sample):
	with open(trim_out_file) as outfile:
		for line in outfile:
			line_split=line.split()
			if line_split[0].startswith("Input"):
				num_reads_input=line_split[3]
				print num_reads_input
				num_reads_surviving=line_split[6]
				print num_reads_surviving
				perc_reads_surviving=line_split[7][1:-2]
				print perc_reads_surviving
				sample_dictionary[sample]=[num_reads_input,num_reads_surviving,perc_reads_surviving]
	return sample_dictionary

def trim_table(sample_dictionary):
	trim_table_filename = "/home/ljcohen/osmotic/evaluation_data/"+"trim_reads_data.txt"
	header=["Sample","Input Reads","Surviving Reads","Percent Surviving"]
    	with open(trim_table_filename,"w") as datafile:
        	datafile.write("\t".join(header))
        	datafile.write("\n")
        	for sample in sample_dictionary.keys():
            		important_nums=sample_dictionary[sample]
            		datafile.write(sample+"\t")
            		datafile.write("\t".join(important_nums))
            		datafile.write("\n")
    	datafile.close()
    	print "Trimmomatic stats written:",trim_table_filename


def execute(listoffiles,basedir):
	sample_dictionary = {}
	pairs_dictionary=get_pairs(listoffiles,basedir)
	#print pairs_dictionary
	for sample in pairs_dictionary:
		print sample
		if sample.startswith("F_heteroclitus"):
			trim_log_dir = "/home/ljcohen/osmotic_trim_F_heteroclitus/"
		else:
			trim_log_dir = "/home/ljcohen/osmotic_trim/trim_log/"
		trim_out_file=trim_log_dir+"trim."+sample+".log"
		matching_string = "TrimmomaticPE: Completed successfully"
        	if os.path.isfile(trim_out_file):
                	with open(trim_out_file) as f:
                        	content = f.readlines()
                	trim_complete = [m for m in content if matching_string in m]
                	if len(trim_complete)!=0:
                        	print "Already trimmed."
				sample_dictionary = get_sample_dictionary(sample_dictionary,trim_out_file,sample)
			else:
				print "Trimmomatic ran, but did not complete:",trim_out_file
		else:
			print "Trim file does not exist:",sample
	trim_table(sample_dictionary)

basedir="/home/ljcohen/osmotic_combined/"
trimdir = "/home/ljcohen/osmotic_trim/"
listoffiles = os.listdir(basedir)
execute(listoffiles,basedir)
