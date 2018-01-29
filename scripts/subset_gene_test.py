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
import time

sys.path.append(os.path.dirname(os.getcwd()))

from models.nn_gtex import MLP
from GTEx import GTEx

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
	top_15 = sort[len(sort) - 15:]
	prev_combos = top_15[:,0].tolist()

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

# get random gene indexes between 0-56238
def create_random_subset(num_genes, tot_gene_lists):		
	#Generate Gene Indexes for Random Sample
	gene_indexes = np.random.randint(0, 56238, num_genes)
	return [tot_gene_lists[i] for i in gene_indexes]


def load_data(num_samples_json, gtex_gct_flt):
	sample_count_dict = {}
	with open(num_samples_json) as f:
		sample_count_dict = json.load(f)

	idx = 0
	data = {}

	for k in sorted(sample_count_dict.keys()):
		data[k] = gtex_gct_flt[:,idx:(idx + sample_count_dict[k])]
		idx = idx + sample_count_dict[k]

	return data


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Run tests on specified subsets of a hallmark or random set')
	parser.add_argument('--set', help='subset to be used', type=str, required=True, choices=['hedgehog', 'notch', 'random'])
	parser.add_argument('--num_genes', help='number of genes', type=int, required=True)
	args = parser.parse_args()

	print('loading genetic data...')
	gtex_gct_flt = np.load('../datasets/gtex_gct_data_float.npy')
	total_gene_list = np.load('../datasets/gtex_complete_gene_list_str.npy')
	print('done')

	data = load_data("../data_scripts/numsamples.json", gtex_gct_flt)

	# load the hedgehog data
	if args.set == 'hedgehog':
		sub = np.load('../datasets/hallmark_numpys/HALLMARK_HEDGEHOG_SIGNALING.npy')
		genes = sub[:,1].tolist()
	elif args.set == 'notch':
		sub = np.load('../datasets/hallmark_numpys/HALLMARK_NOTCH_SIGNALING.npy')
		genes = sub[:,1].tolist()
	else:
		genes = create_random_subset(args.num_genes, total_gene_list)


	print('beginning search for optimal combinations...')
	for i in xrange(1, len(genes)):
		print('--------ITERATION ' + str(i) + '--------')
		# read in the previous accuracy file
		if i > 3:
			# for combos from files
			gene_dict = create_new_combos_from_file('../logs/hedgehog/hh_' + str(i - 1) + '_gene_accuracy.txt', genes)
			# create files to write to, specify neural net architecture
			files = ['hh_' + str(i) + '_gene_accuracy.txt']
		else:
			# for all possible combos
			gene_dict = create_raw_combos(genes, i)
			
			# create files to write to
			files = ['hh_' + str(i) + '_gene_accuracy.txt']
		
		# define hidden layer sizes
		h1 = [1024]
		h2 = [1024]
		h3 = [1024]

		# open log file to write to
		fp = open('../logs/hedgehog/' + files[0], 'w')

		for key in gene_dict:
			# retrieve the new combination of genes and create a new dataset containing the specified features
			start = time.clock()
			combo = list(key)
			#create_subset(combo, total_gene_list)

			gtex = GTEx(data, total_gene_list, combo)

			# partition the newly created datset into a training and test set
			os.system('python ../data_scripts/create-sets.py -d gtex -p ../datasets/GTEx_Rand ' + ' -t 70 -r 30 ')
			
			# run the neural network architecture to retrieve an accuracy based on the new dataset
			mlp = MLP(n_input=i, n_classes=53, batch_size=256, lr=0.001, epochs=75, n_h1=h1[0], n_h2=h2[0], n_h3=h3[0])
			acc = mlp.run(gtex)

			print(str(combo) + '\t' + str(acc))
			
			fp.write('{0}\t{1}\n'.format(key, acc))

		fp.close() 
