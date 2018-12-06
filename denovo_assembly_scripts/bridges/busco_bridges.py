import os
import os.path
import pandas as pd
# custom Lisa module
import clusterfunc

def run_busco(fasta,sample,buscodir):
    busco_command="""
source ~/.bashrc
source activate py3.dammit
export DAMMIT_DB_DIR=$SCRATCH/dammit
cd /pylon5/bi5fpmp/ljcohen/kfish_busco
run_BUSCO.py -i {} -o {} -l /pylon5/bi5fpmp/ljcohen/dammit/busco2db/metazoa_odb9 -m tran --cpu 4
""".format(fasta,sample)
    print(busco_command)
    commands = [busco_command]
    process_name = "busco"
    module_name_list = ""
    filename = sample
    clusterfunc.sbatch_file(buscodir, process_name,module_name_list, filename, commands)

def parse_busco_stats(busco_filename, sample):
    count = 0
    important_lines = [10, 13, 14, 15]
    busco_dict = {}
    busco_dict[sample] = []
    if os.path.isfile(busco_filename):
    	if os.stat(busco_filename).st_size != 0:
        	with open(busco_filename) as buscofile:
            		for line in buscofile:
                                count += 1
                                #print(count)
                                line_data = line.split()
                                #print(line_data)
                                if count in important_lines:
                                    busco_dict[sample].append(int(line_data[0]))
    busco_data = pd.DataFrame.from_dict(busco_dict, orient='index')
    busco_data.columns = ["Complete", "Fragmented", "Missing", "Total"]
    busco_data['Complete_BUSCO_perc'] = busco_data[
        'Complete'] / busco_data['Total']
    return busco_data

def build_DataFrame(data_frame, busco_data):
    # columns=["sample","Complete","Fragmented","Missing","Total"]
    frames = [data_frame, busco_data]
    data_frame = pd.concat(frames)
    return data_frame

def execute(fasta_files,basedir,busco_dir,data_frame):
    count = 0
    # construct an empty pandas dataframe to add on each assembly.csv to
    for filename in fasta_files:
        if filename.endswith(".fasta"):
            sample= filename.split(".")[0]
            print(sample)
            fasta = basedir + filename
            busco_file = busco_dir + "run_" + sample + "/short_summary_" + sample + ".txt"
            if os.path.isfile(busco_file):
                count += 1
                data = parse_busco_stats(busco_file, sample)
                data_frame = build_DataFrame(data_frame, data)
            else:
                run_busco(fasta,sample,busco_dir) 
    return data_frame

basedir = "/pylon5/bi5fpmp/ljcohen/kfish_trinity/"
busco_dir = "/pylon5/bi5fpmp/ljcohen/kfish_busco/"
data_frame = pd.DataFrame()
fasta_files = os.listdir(basedir)
data_frame = execute(fasta_files,basedir,busco_dir,data_frame)
data_frame.to_csv("../../evaluation_data/busco_scores_Dec2018_metazoa.csv")
print("File written: ../../evaluation_data/busco_scores_Dec2018_metazoa.csv")

