import os
import pandas as pd

basedir = "/home/nreid/rnaseq/rawdata/"
flowcell_dirs = os.listdir(basedir)
sample_dict = {}
unique_lanes = []
for flowcell in flowcell_dirs:
    flowcell_dir = basedir + flowcell + "/"
    samples = os.listdir(flowcell_dir)
    for sample in samples:
        if sample.startswith("Sample"):
            sample_info = sample.split("_")
            pool = sample_info[1]
            if sample_info[3] == "heteroclitus":
                genus_species = sample_info[2] + "_" + sample_info[3] + sample_info[4]
                treatment = sample_info[5]
                replicate = sample_info[6]
            else:
                genus_species = sample_info[2]+"_"+sample_info[3]
                treatment = sample_info[4]
                replicate = sample_info[5]
            sample_rep_info = genus_species + "_" + treatment + "_" + replicate
            sample_dir = flowcell_dir + sample + "/"
            files = os.listdir(sample_dir)
            for filename in files:
                if filename.endswith(".gz"):
                    files_info = filename.split("_") 
                    if files_info[2] == "heteroclitus":
                        barcode = files_info[6]
                        lane = files_info[7]
                    else:
                        barcode = files_info[5]
                        lane = files_info[6]
                        read = files_info[7]
                        file_replicate = files_info[8]
                        sample_info_tuple = (flowcell,pool,lane,barcode)
                    if sample_rep_info in sample_dict:
                        if sample_info_tuple not in sample_dict[sample_rep_info]:
                            sample_dict[sample_rep_info].append(sample_info_tuple)
                    else:
                        sample_dict[sample_rep_info]=[sample_info_tuple]
                    pool_lane_info = flowcell + "_" + pool + "_" + lane
                    if pool_lane_info not in unique_lanes:
                        unique_lanes.append(pool_lane_info)
#data_frame = pd.DataFrame({'Sample':'NA','Info':'NA'})
data_frame = pd.DataFrame()
for sample_name in sample_dict:
    print(sample_name)
    for info in sample_dict[sample_name]:
       print(info)
       df2 = pd.DataFrame({'Sample':sample_name,'Flowcell':[info[0]],'Pool':[info[1]],'Lane':[info[2]],'Barcode':[info[3]]})
       data_frame = pd.concat([df2, data_frame])
print(sample_dict.keys())
print("Num samples:",len(sample_dict.keys()))
print(unique_lanes)
print("Num unique lanes:",len(unique_lanes))
data_frame.to_csv("osmotic_killifish_RNAseq_sample_info.csv")
print("File written: osmotic_killifish_RNAseq_sample_info.csv")
