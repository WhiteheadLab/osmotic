import os
import os.path
import pandas as pd
# custom Lisa module
import clusterfunc

def parse_transrate_stats(transrate_assemblies):
    print(transrate_assemblies)
    if os.stat(transrate_assemblies).st_size != 0:
        data = pd.DataFrame.from_csv(transrate_assemblies, header=0, sep=',')
        return data

def build_DataFrame(data_frame, transrate_data):
    # columns=["n_bases","gc","gc_skew","mean_orf_percent"]
    frames = [data_frame, transrate_data]
    data_frame = pd.concat(frames)
    return data_frame

def transrate_forward(transratedir,transrate_out,trinity_fasta,sample,reference):
    transrate_command = """
source ~/.bashrc
transrate --assembly={} --reference={} \
--threads=8 \
--output={}
""".format(trinity_fasta,reference,transrate_out)
    commands = [transrate_command]
    process_name = "transrate_forward"
    module_name_list = ""
    filename = sample + "_forward"
    clusterfunc.sbatch_file(transratedir, process_name,module_name_list, filename, commands)

def transrate_reverse(transratedir,transrate_out,trinity_fasta,sample,reference):
    transrate_command = """
source ~/.bashrc
transrate --assembly={} --reference={} \
--threads=8 \
--output={}
""".format(reference,trinity_fasta,transrate_out)
    commands = [transrate_command]
    process_name = "transrate_reverse"
    module_name_list = ""
    filename = sample + "_reverse"
    clusterfunc.sbatch_file(transratedir, process_name,module_name_list, filename, commands)

def execute(data_frame1,data_frame2,new_assemblies,old_assemblies,new_assemblies_dir,old_assemblies_dir,transrate_dir,reference):
    # construct an empty pandas dataframe to add on each assembly.csv to
    for fasta in new_assemblies:
        if fasta.endswith(".fasta"):
            print(fasta)
            sample = fasta.split(".")[0]
            trinity_fasta = new_assemblies_dir + fasta 
            transrate_out_forward = transrate_dir + sample + "_trinity_v_Fhet.NCBIrna/"
            transrate_out_reverse = transrate_dir + sample + "_Fhet.NCBIrna_v_trinity/"
            transrate_assemblies_forward = transrate_out_forward + "assemblies.csv"
            transrate_assemblies_reverse = transrate_out_reverse + "assemblies.csv"
            if os.path.isfile(transrate_assemblies_forward):
                data1 = parse_transrate_stats(transrate_assemblies_forward)
                data_frame1 = build_DataFrame(data_frame1, data1)
            else:  
                print("Forward transrate output not found:",transrate_assemblies_forward)
                print("Running transrate forward...")
                transrate_forward(transrate_dir,transrate_out_forward,trinity_fasta,sample,reference)
            if os.path.isfile(transrate_assemblies_reverse):
                data2 = parse_transrate_stats(transrate_assemblies_reverse)
                data_frame2 = build_DataFrame(data_frame2,data2)
            else:
                print("Reverse transrate output not found:",transrate_assemblies_reverse)
                print("Running transrate reverse...")
                transrate_reverse(transrate_dir,transrate_out_reverse,trinity_fasta,sample,reference)
    return data_frame1, data_frame2


reference = "/pylon5/bi5fpmp/ljcohen/Fhet_reference_genome/ncbi/rna.fa"
#reference = "/pylon5/bi5fpmp/ljcohen/Fhet_reference_genome/ensembl/Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.cds.all.fa"
#reference = "/pylon5/bi5fpmp/ljcohen/Fhet_reference_genome/ensembl/Fundulus_heteroclitus.Fundulus_heteroclitus-3.0.2.cdna.all.fa"
transrate_dir = "/pylon5/bi5fpmp/ljcohen/kfish_transrate/"
old_assemblies_dir = "/pylon5/bi5fpmp/ljcohen/kfish_assemblies_old/"
new_assemblies_dir = "/pylon5/bi5fpmp/ljcohen/kfish_trinity/"
old_assemblies = os.listdir(old_assemblies_dir)
new_assemblies = os.listdir(new_assemblies_dir)
data_frame1 = pd.DataFrame()
data_frame2 = pd.DataFrame()
data_frame1,data_frame2 = execute(data_frame1,data_frame2,new_assemblies,old_assemblies,new_assemblies_dir,old_assemblies_dir,transrate_dir,reference)
data_frame1.to_csv("../../evaluation_data/transrate_reference_trinity_v_Fhet.NCBI_Dec2018.csv")
data_frame2.to_csv("../../evaluation_data/transrate_reference_Fhet.NCBI_v_trinity_Dec2018.csv")
#print("Forward transrate reference written: evaluation_data/transrate_reference_trinity_v_Fhet.NCBI.csv")
#print("Reverse transrate reference written: evaluation_data/transrate_reference_Fhet.NCBI_v_trinity.csv")
