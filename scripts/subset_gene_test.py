#/usr/bin/python

'''
	This file takes divides a hallmark subset into user-specified subsets in order to
	determine the best combination of genes for classification purposes.
'''

import numpy as np 
import os
import subprocess
import json
import shutil

sub = np.load('../datasets/hallmark_numpys/HALLMARK_HEDGEHOG_SIGNALING.npy')

genes = sub[:,1].tolist()

data = sub[:,2:] # raw data is in 2

def create_subset(sub_gene, tot_gene_lists):
	gene_count_dict = dict()
	with open("../data_scripts/numsamples_gtex.json") as f:
	    gene_count_dict = json.load(f)

	gene_indexes = []
	gene_indexes.append(np.argwhere(tot_gene_lists == sub_gene))

	if os.path.isdir('../datasets/GTEx_Hedgehog'):
		shutil.rmtree('../datasets/GTEx_Hedgehog')

	if not os.path.isdir('./datasets/GTEx_Hedgehog'):
		os.mkdir('../datasets/GTEx_Hedgehog')

	tissues = os.listdir("../datasets/GTEx_Data")


	for i in range(len(tissues)):#loop through all tissue types
	    tissue_samples = os.listdir("../datasets/GTEx_Data/"+ tissues[i])#get dat files/samples

	    if not os.path.isdir('./datasets/GTEx_Hedgehog/' + tissues[i]):
	    	os.mkdir('../datasets/GTEx_Hedgehog/' + tissues[i])

	    for j in range(gene_count_dict[tissues[i]]):#loop through dat files
	        sample = np.fromfile("../datasets/GTEx_Data/"+ tissues[i]+'/' + tissue_samples[j], dtype=np.float32)#get 1 sample
		random_sample = np.zeros(1)
		for r in range(len(gene_indexes)):#loop through random indexes
			random_sample[r] = sample[gene_indexes[r]]

		random_sample = random_sample.astype(np.float32)
		#print(len(random_sample))
		random_sample.tofile("../datasets/GTEx_Hedgehog/"+ tissues[i] + '/' + tissue_samples[j])
		



files = ['hedgehog_single_gene_accuracy.txt']
h1 = [1024]
h2 = [1024]
h3 = [1024]

accs = np.zeros((len(genes),1))

gtex_str_data = np.load('../datasets/gtex_gct_data_string.npy')
total_gene_list = gtex_str_data[1:,1]

j = 0

for i in range(len(genes)):
	#os.system('python ../data_scripts/random_gene_generator.py --num ' + str(1))
	create_subset(genes[i], total_gene_list)

	os.system('python ../data_scripts/create-sets.py -d gtex -p ../datasets/GTEx_Hedgehog ' + ' -t 70 -r 30 ')
	#find features

	#future:change arguments and get from command line
	accs[i,j] = subprocess.check_output('python ../models/nn_gtex.py --n_input ' + str(1) + \
		' --n_classes 53 --batch_size 256 --lr 0.001 --epochs 75 --h1 ' + str(h1[0]) + ' --h2 ' + str(h2[0]) + ' --h3 ' + str(h3[0]), shell=True)

np.save('../logs/accuracies_hedgehog.npy', accs)

means = np.mean(accs, axis=1)

fp = open('../logs/' + files[0], 'w')
for val,name in zip(means,genes):
    fp.write('{0}\t{1}\n'.format(name, val))    
fp.close() 
