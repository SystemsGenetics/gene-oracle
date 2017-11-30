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
import sys, argparse

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
def create_new_combos_from_file(file, genes):
	combos = []
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

# create every possible combination
def create_raw_combos(genes, i):
	combos = []
	for c in itertools.combinations(genes, i):
		combos.append(c)

	return dict.fromkeys(combos)

# create a new dataset in the datasets directory that contains the given genes in sub_gene
def create_subset(sub_gene, tot_gene_lists):
	gene_count_dict = dict()
	with open("../data_scripts/numsamples_gtex.json") as f:
	    gene_count_dict = json.load(f)

	gene_indexes = []
	for i in range(len(sub_gene)):
		gene_indexes.append(np.argwhere(tot_gene_lists == sub_gene[i])[0])

	if os.path.isdir('../datasets/GTEx_Notch'):
		shutil.rmtree('../datasets/GTEx_Notch')

	if not os.path.isdir('./datasets/GTEx_Notch'):
		os.mkdir('../datasets/GTEx_Notch')

	tissues = os.listdir("../datasets/GTEx_Data")


	for i in range(len(tissues)):#loop through all tissue types
	    tissue_samples = os.listdir("../datasets/GTEx_Data/"+ tissues[i])#get dat files/samples

	    if not os.path.isdir('./datasets/GTEx_Notch/' + tissues[i]):
	    	os.mkdir('../datasets/GTEx_Notch/' + tissues[i])

	    for j in range(gene_count_dict[tissues[i]]):#loop through dat files
	        sample = np.fromfile("../datasets/GTEx_Data/"+ tissues[i]+'/' + tissue_samples[j], dtype=np.float32)#get 1 sample

		sub_sample = np.zeros(len(gene_indexes))
		for r in range(len(gene_indexes)):#loop through indexes
			sub_sample[r] = sample[gene_indexes[r]]

		sub_sample = sub_sample.astype(np.float32)
		#print(len(random_sample))
		sub_sample.tofile("../datasets/GTEx_Notch/"+ tissues[i] + '/' + tissue_samples[j])

def create_random_subset(num_genes, tot_gene_lists):		
	#Generate Gene Indexes for Random Sample
	gene_indexes = np.random.randint(0, 56238, num_genes)
	return [tot_gene_lists[i] for i in gene_indexes]


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Run tests on specified subsets of a hallmark or random set')
	parser.add_argument('--set', help='subset to be used', type=str, required=True, choices=['hedgehog', 'notch', 'random'])
	parser.add_argument('--num_genes', help='number of genes', type=int, required=True)
	args = parser.parse_args()

	gtex_str_data = np.load('../datasets/gtex_gct_data_string.npy')
	total_gene_list = gtex_str_data[1:,1]

	# load the hedgehog data
	if args.set == 'hedgehog':
		sub = np.load('../datasets/hallmark_numpys/HALLMARK_HEDGEHOG_SIGNALING.npy')
		genes = sub[:,1].tolist()
	elif: args.set == 'notch':
		sub = np.load('../datasets/hallmark_numpys/HALLMARK_NOTCH_SIGNALING.npy')
		genes = sub[:,1].tolist()
	else:
		genes = create_random_subset(args.num_genes, total_gene_list)

	for i in xrange(1, len(genes)):

		print('--------ITERATION ' + str(i) + '--------')
		# read in the previous accuracy file
		if i > 3:
			gene_dict = create_new_combos_from_file('../logs/notch/notch_' + str(i) + '_gene_accuracy.txt', genes)
		else:
			gene_dict = create_raw_combos(genes, i)

		# collect the data from the hedgehog file
		data = sub[:,2:] # raw data is in 2

		# create files to write to, specify neural net architecture
		files = ['notch_' + str(i + 1) + '_gene_accuracy.txt']
		
		h1 = [1024]
		h2 = [1024]
		h3 = [1024]

		fp = open('../logs/notch/' + files[0], 'w')

		for key in gene_dict:
			# retrieve the new combination of genes and create a new dataset containing the specified features
			combo = list(key)
			create_subset(combo, total_gene_list)

			# partition the newly created datset into a training and test set
			os.system('python ../data_scripts/create-sets.py -d gtex -p ../datasets/GTEx_Notch ' + ' -t 70 -r 30 ')

			# run the neural network architecture to retrieve an accuracy based on the new dataset
			acc = subprocess.check_output('python ../models/nn_gtex.py --n_input ' + str(i + 1) + \
				' --n_classes 53 --batch_size 256 --lr 0.001 --epochs 75 --h1 ' + str(h1[0]) + ' --h2 ' + str(h2[0]) + ' --h3 ' + str(h3[0]), shell=True)

			print('iteration ' + str(i) + ' ' + str(combo) + '\t' + str(acc))
			
			fp.write('{0}\t{1}\n'.format(key, acc))

		fp.close() 
