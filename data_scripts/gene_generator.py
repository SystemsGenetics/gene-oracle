import numpy as np
import os
import json
import shutil
import sys, argparse

parser = argparse.ArgumentParser(description='generate random subsets of GTEx specified by number')
parser.add_argument('--num', help='number of genes to use in random submatrix', type=int, required=True, default=0)
args = parser.parse_args()


#Load JSON
gene_count_dict = dict()
with open("numsamples_gtex.json") as f:
    gene_count_dict = json.load(f)

if os.path.isdir('../datasets/GTEx_Data_Random'):
	shutil.rmtree('../datasets/GTEx_Data_Random')

if not os.path.isdir('./datasets/GTEx_Data_Random'):
	os.mkdir('../datasets/GTEx_Data_Random')

tissues = os.listdir("../datasets/GTEx_Data")

for i in range(len(tissues)):#loop through all tissue types
    tissue_samples = os.listdir("../datasets/GTEx_Data/"+ tissues[i])#get dat files/samples

    if not os.path.isdir('./datasets/GTEx_Data_Random/' + tissues[i]):
    	os.mkdir('../datasets/GTEx_Data_Random/' + tissues[i])

    for j in range(gene_count_dict[tissues[i]]):#loop through dat files
        sample = np.fromfile("../datasets/GTEx_Data/"+ tissues[i]+'/' + tissue_samples[j], dtype=np.float32)#get 1 sample
	random_sample = np.zeros(args.num)
	for r in range(len(gene_indexes)):#loop through random indexes
		random_sample[r] = sample[gene_indexes[r]]

	random_sample = random_sample.astype(np.float32)
	#print(len(random_sample))
	random_sample.tofile("../datasets/GTEx_Data_Random/"+ tissues[i] + '/' +tissue_samples[j])
		
        
