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
import itertools

# special helper function to sanitize the string containing the genes from
# an accuracy file
def sanitize(gene_str):
	gene_list = gene_str.strip('()')
	gene_list = gene_list.replace('\'', '')
	gene_list = gene_list.replace(' ', '')
	gene_list = gene_list.split(',')
	return gene_list

# create_new_combos takes in a file string that is an accuracy file with a list of genes
# separated by a tab, followed by the accuracy for that list. it returns a dictionary
# of new combinations with one extra gene appended that was not previously in the list
def create_new_combos(file, genes):
	combos = []
	# for subset in itertools.combinations(genes, 2):
	# 	combos.append(subset)
	prev_accs = np.loadtxt(file, delimiter='\t', dtype=np.str)
	sort = prev_accs[np.argsort(prev_accs[:,1])]
	top_5_perc = sort[len(sort) - 15:]
	prev_combos = top_5_perc[:,0].tolist()

	for c in prev_combos:
		gene_list = sanitize(c)

		for g in genes:
			if g not in gene_list:
				temp_list = gene_list[:]
				temp_list.append(g)
				combos.append(temp_list)

	combos = [tuple(g) for g in combos]
	return dict.fromkeys(combos)

# create a new dataset in the datasets directory that contains the given genes in sub_gene
def create_subset(sub_gene, tot_gene_lists):
	gene_count_dict = dict()
	with open("../data_scripts/numsamples_gtex.json") as f:
	    gene_count_dict = json.load(f)

	gene_indexes = []
	for i in range(len(sub_gene)):
		gene_indexes.append(np.argwhere(tot_gene_lists == sub_gene[i]))

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

		sub_sample = np.zeros(len(gene_indexes))
		for r in range(len(gene_indexes)):#loop through random indexes
			sub_sample[r] = sample[gene_indexes[r]]

		sub_sample = sub_sample.astype(np.float32)
		#print(len(random_sample))
		sub_sample.tofile("../datasets/GTEx_Hedgehog/"+ tissues[i] + '/' + tissue_samples[j])
		

# load the hedgehog data
sub = np.load('../datasets/hallmark_numpys/HALLMARK_HEDGEHOG_SIGNALING.npy')
genes = sub[:,1].tolist()

for i in xrange(6, 36):

	# read in the previous accuracy file
	gene_dict = create_new_combos('../logs/hedgehog_' + str(i) + '_gene_accuracy.txt', genes)

	# collect the data from the hedgehog file
	data = sub[:,2:] # raw data is in 2

	# create files to write to, specify neural net architecture
	files = ['hedgehog_' + str(i + 1) + '_gene_accuracy.txt']
	h1 = [1024]
	h2 = [1024]
	h3 = [1024]

	# collect the total list of genes in the GTEx dataset
	gtex_str_data = np.load('../datasets/gtex_gct_data_string.npy')
	total_gene_list = gtex_str_data[1:,1]

	fp = open('../logs/' + files[0], 'w')
	i = 0
	for key in gene_dict:
		# retrieve the new combination of genes and create a new dataset containing the specified features
		combo = list(key)
		create_subset(combo, total_gene_list)

		# partition the newly created datset into a training and test set
		os.system('python ../data_scripts/create-sets.py -d gtex -p ../datasets/GTEx_Hedgehog ' + ' -t 70 -r 30 ')

		# run the neural network architecture to retrieve an accuracy based on the new dataset
		acc = subprocess.check_output('python ../models/nn_gtex.py --n_input ' + str(6) + \
			' --n_classes 53 --batch_size 256 --lr 0.001 --epochs 75 --h1 ' + str(h1[0]) + ' --h2 ' + str(h2[0]) + ' --h3 ' + str(h3[0]), shell=True)

		print('iteration ' + str(i) + ' ' + str(combo) + '\t' + str(acc))
		
		fp.write('{0}\t{1}\n'.format(key, acc))

		i = i + 1

	fp.close() 
