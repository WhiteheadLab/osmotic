import pandas as pd
import os
from itertools import izip

def get_annotation_counts(species_annotation_file,counts_data):
    annotation = pd.read_table(species_annotation_file,header=0)
    annotation = annotation.drop(["Start","End","Strand"],axis=1)
    frames = [annotation,counts_data]
    annotated_counts = pd.concat(frames,axis=1,join='inner')
    annotated_counts = annotated_counts.drop(["Chr"],axis=1)
    return annotated_counts

# problem is there are more than 1 transcript ID per GeneID
# Check for uniqueness
# if more than 1:
# if there 0 counts?
# drop 0 counts
# elif more than one transcript ID counted per GeneID?
# then get sum of different transcript ID, same GeneID

def get_unique_annotations(annotated_counts_old):
    # assume rows are sorted, with GeneID and transcriptID adjacent
    # is this a safe assumption?
    # probably not
    # sort
    annotated_itterows = annotated_counts_old.iterrows()
    check_GeneID = ""
    check_transcript = ""
    for i, row in annotated_itterows:
        GeneID = (row['GeneID'])
        transcript = (row['transcript'])
        if GeneID == check_GeneID:
            if row['count'] == 0:
                annotated_counts_old.drop(i)
            elif row['count'] > 0:
                # add counts together
                new_count += count
                # append count value in current row
                annotated_counts_old.set_value(i[1],row['count'], new_count)
                # drop previous row
                annotated_counts_old.drop(previous_i)
        else:
            new_count = 0
        check_GeneID = GeneID
        check_transcript = transcript
        count = row['count']
        previous_row = row
        previous_i = i
    return annotated_counts_old
            
def parse_data(counts_file):
    counts_data=pd.read_table(counts_file,header=0)
    return counts_data
    
def build_DataFrame(counts_dataframe,counts_data):
    counts_dataframe.head()
    counts_data.head()
    counts_dataframe=pd.merge(counts_dataframe,counts_data,on=["GeneID"])
    return counts_dataframe


def execute(dirnames,annotation_dir,salmon_dir):
	with open("/home/ljcohen/osmotic/kfish2rae5g.GeneID","rU") as GeneID_file:
    		GeneID = pd.read_table(GeneID_file,header=0)
	counts_dataframe = GeneID
	for genus_species in dirnames:
        	counts_dir = salmon_dir+genus_species
        	counts_files = os.listdir(counts_dir)
        	species_annotation_dir = annotation_dir+genus_species+"/"
        	species_annotation_file = species_annotation_dir+genus_species+".saf"
        	print species_annotation_file 
        	for filename in counts_files:
            		if filename.endswith(".counts"):
               			new_filename=filename.split(".")[0]
                		print new_filename
				#counts_file=counts_dir+"/"+filename
                		#counts_data=parse_data(counts_file)
                		#annotated_counts = get_annotation_counts(species_annotation_file,counts_data)
                		#old_annotation_counts = get_unique_annotations(annotated_counts)
                		#counts_dataframe=build_DataFrame(counts_dataframe,annotated_counts)
                		#counts_dataframe.rename(columns = {'count': new_filename},inplace=True)
        

salmon_dir="/home/ljcohen/osmotic_salmon/"
annotation_dir="/home/ljcohen/msu_assemblies_finished/"
dirnames=os.listdir(salmon_dir)
print(dirnames)
execute(dirnames,annotation_dir,salmon_dir)

