import numpy as np
import os
import json

#Generate Gene Indexes for Random Sample
gene_indexes = np.random.randint(56238, size=200)

#Load JSON
gene_count_dict = dict()
with open("numsamples_gtex.json") as f:
    gene_count_dict = json.load(f)

tissues = os.listdir("./datasets/GTEx_Data")

for i in range(len(tissues)):#loop through all tissue types
    tissue_samples = os.listdir("./datasets/GTEx_Data/"+ tissues[i])#get dat files/samples

    for j in range(gene_count_dict[tissues[i]]):#loop through dat files
        sample = np.fromfile("./datasets/GTEx_Data/"+ tissues[i]+'/' + tissue_samples[j], dtype=np.float32)#get 1 sample
	random_sample = np.zeros(200)
	for r in range(len(gene_indexes)):#loop through random indexes
		random_sample[r] = sample[gene_indexes[r]]

	random_sample = random_sample.astype(np.float32)
	#print(len(random_sample))
	random_sample.tofile("./datasets/GTEx_Data_random200/"+ tissues[i] + '/' +tissue_samples[j])
		
        
