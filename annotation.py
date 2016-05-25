import os
from os import path
import subprocess
from subprocess import Popen, PIPE
import clusterfunc


def create_bed_file(species_dir,species):
	make_bed="""
BT=/home/nreid/bin/bedtools/bin/bedtools2
$BT bamtobed -i {}{}.bam | sortBed -i - > {}{}.sorted.bed
""".format(species_dir,species,species_dir,species)
	return make_bed


def create_id_file(species_dir,species,sorted_bed):
	id_command="""
BT=/home/nreid/bin/bedtools/bin/bedtools2
$BT map -a <(sort -V -k 1,1 -k 2,2n {}) \\
-b <(awk '$3 ~ /CDS/' /home/nreid/popgen/kfish3/kfish2rae5g.main.pub.sort.gff | grep -v scaffold_1167) \\
-c 9 -o distinct,count \\
-g /home/nreid/popgen/kfish3/bed/killifish20130322asm.genome | \\
grep Parent | \\
awk '$7 !~ /,/' | \\
sed 's/Parent=//' | \\
sed 's/;[^\t]*//' > {}{}.id
""".format(sorted_bed,species_dir,species)
	return id_command

def get_transcript_id(id_file):
	transcript_id_data={}
	with open(id_file,"rU") as idfile:
		for line in idfile:
			line_info=line.split("\t")
			transcript_id=line_info[3]
			ref_id=line_info[6]
			if transcript_id not in transcript_id_data.keys():
				transcript_id_data[transcript_id]=ref_id
				#if ref_id not in transcript_id_data[transcript_id]:
				#	transcript_id_data[transcript_id].append(ref_id)
				#	print "Multiple reference ID found:",transcript_id,transcript_id_data[transcript_id]
			#else:
			#	transcript_id_data[transcript_id]=[ref_id]
	print "Number of unique contigs:",len(transcript_id_data.keys())
	return transcript_id_data	

basedir="/home/ljcohen/msu_assemblies_finished/"
#listofdirs=os.listdir(basedir)
listofdirs=["F_diaphanus/F_diaphanus.trinity.2","F_sciadicus/F_sciadicus.trinity.2"]
trinityfile="Trinity.fasta"
for dirs in listofdirs:
	#if dirname == "F_chrysotus":
		dirname = dirs.split("/")[0]
		newdir=basedir+dirs+"/"
		trinity=newdir+trinityfile
		print trinity
		id_file=newdir+dirname+".id"
		sorted_bed_file=newdir+dirname+".sorted.bed"
		print sorted_bed_file
		make_bed=create_bed_file(newdir,dirname)
		print make_bed
		#s=subprocess.Popen(make_bed,shell=True)
		#s.wait()
		#print id_file
		id_command=create_id_file(newdir,dirname,sorted_bed_file)
		print id_command
		id_command=[id_command]
		process_name="annotation_id"
		module_name_list=""
		#clusterfunc.sbatch_file(newdir,process_name,module_name_list,dirname,id_command)
		saf_file=newdir+dirname+".saf"
		transcript_id_data=get_transcript_id(id_file)
		with open(trinity,"rU") as transcriptome_assembly:
			sequence=""
			with open(saf_file,"w") as saf:
				saf.write("GeneID"+"\t")
				saf.write("Chr"+"\t")
				saf.write("Start"+"\t")
				saf.write("End"+"\t")
				saf.write("Strand"+"\n")
				for line in transcriptome_assembly:
					if line.startswith(">"):
						sequence=""
						split_line=line.split(" ")
						transcript_id=split_line[0][1:]
					else:
						if len(line) == 61: 
							sequence+=line
						elif len(line) < 61 and len(line) > 0: 
							sequence+=line
							contig_length=len(sequence)
							print contig_length
							if transcript_id in transcript_id_data.keys():
								ref_id=transcript_id_data[transcript_id]
								print transcript_id
								print ref_id
								saf.write(ref_id+"\t")
                	                               		saf.write(transcript_id+"\t")
                                                		saf.write("1"+"\t")
								saf.write(str(contig_length)+"\t")
								saf.write("-"+"\n")
						elif len(line) == 0:
							print "Blank line?"
	#else:
	#	print "Skipping:", dirname
