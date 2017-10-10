import os
import json

tissues = os.listdir("./datasets/GTEx_Data")
gene_count_dict = dict()

#Create dictionary of the number of samples in each Tissue Directory
for i in range(len(tissues)):#loop through all tissue types
    tissue_samples = os.listdir("./datasets/GTEx_Data/"+ tissues[i])#get dat files/samples
    for j in range(len(tissue_samples)):
        gene_count_dict[tissues[i]] = len(tissue_samples)

with open("numsamples_gtex.json", 'w') as f:
    json.dump(gene_count_dict, f)

print(sum(gene_count_dict.values()))
