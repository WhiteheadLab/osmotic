import os
import os.path
from os.path import basename
from urllib import urlopen
from urlparse import urlparse
import subprocess
from subprocess import Popen, PIPE
import urllib
import shutil
import glob
# custom Lisa module
import clusterfunc
# Python plotting libraries
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats, integrate

def get_assemblies(assemblydir):
	genus_species_dirs=os.listdir(assemblydir)
	data_frame=pd.DataFrame()
	for genus_species in genus_species_dirs:
		genus_species_dir=assemblydir+genus_species+"/"
		transrate_dir=genus_species_dir+"kf2evg367mixx11.transrate/"
		transrate_assemblies=transrate_dir+"assemblies.csv"
		if os.path.isfile(transrate_assemblies):
			data=parse_transrate_stats(transrate_assemblies)
			data_frame=build_DataFrame(data_frame,data)
		else:
			print "Assemblies not complete:",transrate_assemblies
	return data_frame

def parse_transrate_stats(transrate_assemblies):
        data=pd.DataFrame.from_csv(transrate_assemblies,header=0,sep=',')
        return data

def build_DataFrame(data_frame,transrate_data):
	frames=[data_frame,transrate_data]
	data_frame=pd.concat(frames)
	return data_frame

assemblydir="/home/ljcohen/msu_assemblies_finished/"
data_frame=get_assemblies(assemblydir)
data_frame.to_csv("/home/ljcohen/osmotic/osmotic_killifish_transrate2_data.csv")
print "File written: osmotic_killifish_transrate2_data.csv"
