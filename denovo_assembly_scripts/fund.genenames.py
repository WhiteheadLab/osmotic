import pandas as pd
import os 
# requires dammit env:
# source activate py3.dammit

from dammit.fileio.gff3 import GFF3Parser

dammit_dir = '/home/ljcohen/osmotic_damit/'
dammit_dirs = os.listdir(dammit_dir)
print(dammit_dirs)
for dammit_dirname in dammit_dirs:
    if dammit_dirname != "sbatch_files":
        genus_species = dammit_dirname.split(".")[0]
        dammit_gff = dammit_dir + dammit_dirname + "/" + genus_species + ".trinity_out.Trinity.fasta.dammit.gff3"
        print(dammit_gff)
        annotations = GFF3Parser(filename=dammit_gff).read()
        all_names = annotations.sort_values(by=['seqid'],ascending=True)[['seqid','Name']]
        annotations = annotations.dropna(subset=['Name'])
        fund = annotations[annotations['Name'].str.startswith("gi")]
        names = fund.sort_values(by=['seqid'], ascending=True)[['seqid', 'Name']]
        names_out = '/home/ljcohen/osmotic_assemblies_farm/'+genus_species+'.trinity_out.Trinity.fasta.Fundulus.genenames.csv'
        #names.to_csv(names_out)
        #print("Written:",names_out)
        all_names_out = '/home/ljcohen/osmotic_assemblies_farm/'+genus_species+'.trinity_out.Trinity.fasta.all_gene_names.csv'
        all_names.to_csv(all_names_out)
        print("Written:",all_names_out)
