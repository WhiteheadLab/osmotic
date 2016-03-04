from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio import Entrez

Entrez.email="ljcohen@ucdavis.edu"
handle=Entrez.egquery(db="nuccore",term="Fundulus heteroclitus mitochondrion")
record=Entrez.read(handle)
print record
#gi_list=record["IdList"]
#print gi_list

#record = SeqIO.read("NC_005816.gb","genbank")
#print record
