import os
import scipy.io as sio
import pandas as pd
import urllib

#download Data from GEO
geo_url = "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE110nnn/GSE110823/suppl/GSE110823_RAW.tar"
local_file = "GSE110823_raw.tar"
urllib.request.urlretrieve(geo_url, local_file)
os.system(command = "tar xf {}".format(local_file))
os.system(command = "gzip -d GSM3017261_150000_CNS_nuclei.mat.gz")
local_matrix = os.path.abspath("GSM3017261_150000_CNS_nuclei.mat")

#load data into environment
data = sio.loadmat(local_matrix)

#Digital Expression Matrix
DGE = data['DGE']

#Genes
genes = pd.Series(data['genes']).str.strip(' ')

#Sample types
sample_type = pd.Series(data['sample_type']).str.strip(' ')

#Main cluster assignment
cluster_assignment = pd.Series(data['cluster_assignment']).str.strip(' ')

#Spinal cluster assignment
spinal_cluster_assignment = pd.Series(data['spinal_cluster_assignment']).str.strip(' ')
