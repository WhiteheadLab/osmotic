import os
import os.path
import clusterfunc
from khmer import HLLCounter

def get_unique_kmers(mmetsp,fasta):
    print(fasta)
    counter = HLLCounter(0.1,25)
    counter.consume_seqfile(fasta)
    unique_kmers = counter.estimate_cardinality()
    print(unique_kmers)
    return unique_kmers

def make_unique_kmers_dictionary(sample_dictionary,fasta,species):
    unique_kmers = get_unique_kmers(species,fasta)
    sample_dictionary[species]=unique_kmers
    return sample_dictionary

def make_unique_kmer_table(sample_dictionary,unique_kmers_filename):
    header=["Sample","Unique_kmers"]
    with open(unique_kmers_filename,"w") as datafile:
        datafile.write("\t".join(header))
        datafile.write("\n")
        for sample in sample_dictionary:
            datafile.write(sample+"\t")
            unique_kmers = str(sample_dictionary[sample])
            datafile.write(unique_kmers)
            datafile.write("\n")
    datafile.close()

def execute(fasta_list,unique_kmers_filename,basedir):
    sample_dictionary = {}
    for fasta in fasta_list:
        if fasta.endswith(".fasta"):
            species= fasta.split(".")[0]
            fasta_file = basedir + fasta
            sample_dictionary = make_unique_kmers_dictionary(sample_dictionary,fasta_file,species)
            make_unique_kmer_table(sample_dictionary,unique_kmers_filename)

#basedir = "/pylon5/bi5fpmp/ljcohen/kfish_trinity/"
basedir = "/pylon5/bi5fpmp/ljcohen/kfish_assemblies_old/"
unique_kmers_filename = "/pylon5/bi5fpmp/ljcohen/osmotic/evaluation_data/unique_kmers_kfish_old.txt"
fasta_list = os.listdir(basedir)
execute(fasta_list,unique_kmers_filename,basedir)
print("File written:",unique_kmers_filename)
