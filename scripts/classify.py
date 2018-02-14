#/usr/bin/python

'''
	This script can be used to run a specified dataset, a specified subset of genes,
	or a specified number of random genes for classification. 

	It is required to have a numpy array containing column-wise samples and rows of
	genes. Additionally, it is required to have a numpy vector of genes that are contained
	in the dataset (in the same exact order). 
'''

import numpy as np 
import sys, argparse
import os
import re

sys.path.append(os.path.dirname(os.getcwd()))

from subset_gene_test import create_random_subset, load_data
from models.nn_gtex import MLP
from GTEx import GTEx

# check the arguments are correct for the program
def check_args(args):
	# check dataset is of correct type
	if os.path.exists(args.dataset):
		split = args.dataset.split('.')
		if split[-1] != 'npy':
			print('Dataset file must be a numpy file.')
			sys.exit(1)
	else:
		print('File does not exist!')
		sys.exit(1)

	# check gene list is of correct type
	if os.path.exists(args.gene_list):
		split = args.gene_list.split('.')
		if split[-1] != 'npy':
			print('Gene list file must be a numpy file.')
			sys.exit(1)
	else:
		print('File does not exist!')
		sys.exit(1)

	# check gene list is of correct type
	if os.path.exists(args.sample_json):
		split = args.sample_json.split('.')
		if split[-1] != 'json':
			print('sample file must be a json file.')
			sys.exit(1)
	else:
		print('File does not exist!')
		sys.exit(1)	


# read a csv or txt file that contains a name of a subset followed by a list of genes
def read_subset_file(file):
	with open(file, 'r') as f:
		content = f.readlines()

	# eliminate new line characters
	content = [x.strip() for x in content]

	# split on tabs or commas to create a sublist of set names and genes
	content = [re.split('\t|,', x) for x in content]

	# create a dictionary with keys subset names and values list of genes
	subsets = {}
	for c in content:
		subsets[c[0]] = c[1:]

	return subsets


# perform random classification based on the given parameters
def random_classification(data, total_gene_list, num_genes, iters, out_file):

	accs = []
	f = open(out_file, 'w')
	f.write('Num\tAverage\tStd Dev\n')

	for num in num_genes:
		# set up the neural network
		mlp = MLP(n_input=num, n_classes=len(data), batch_size=128, lr=0.001, epochs=75, n_h1=1024, n_h2=1024, n_h3=1024)
		
		for i in xrange(iters):
			# generate random set of genes from the total gene list
			r_genes = create_random_subset(num_genes, total_gene_list)

			# set up the gtex class to partition data
			gtex = GTEx(data, total_gene_list, r_genes)

			# run the neural net
			accs.append(mlp.run(gtex))

		# print the results to a file
		accs_np = np.asarray(accs)
		mean = np.mean(accs_np)
		std = np.std(accs_np)
		f.write(str(num) + '\t' + str(mean) + '\t' + str(std) + '\n')

	f.close()


# perform classificaiton on each of the subsets provided in the subset_list argument
def subset_classification(data, total_gene_list, subsets, out_file, kfold_val=1):
	accs = []
	f = open(out_file, 'w')
	f.write('Num\tAverage\tStd Dev\n')

	for s in subsets:
		print ('running experiment on ' + str(s))
		# set up the neural network
		mlp = MLP(n_input=len(subsets[s]), n_classes=len(data), batch_size=128, lr=0.001, epochs=75, n_h1=1024, n_h2=1024, n_h3=1024)
		print(s)
		print(subsets[s])		
		for i in xrange(kfold_val):
			# set up the gtex class to partition data
			gtex = GTEx(data, total_gene_list, subsets[s])

			# run the neural net
			accs.append(mlp.run(gtex))

		# print the results to a file
		accs_np = np.asarray(accs)
		mean = np.mean(accs_np)
		std = np.std(accs_np)
		f.write(str(num) + '\t' + str(mean) + '\t' + str(std) + '\n')

	f.close()


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Run classification on specified dataset, \
		subset of genes, or a random set')
	parser.add_argument('--dataset', help='dataset to be used', type=str, required=True)
	parser.add_argument('--gene_list', help='list of genes in dataset (same order as dataset)', \
		type=str, required=True)
	parser.add_argument('--sample_json', help='json file containing number of samples per class', \
		type=str, required=True)
	parser.add_argument('--out_file', help='output file to send results to', type=str, required=True)
	parser.add_argument('--subset_list', help='gmt/gct file containing subsets', type=str, required=False)
	parser.add_argument('--random_test', help='Perform random test', action='store_true', required=False)
	parser.add_argument('--num_random_genes', help='Number of random genes to assess', nargs='+', \
		type=int, required=False)
	parser.add_argument('--rand_iters', help='Number of iterations to perform for random classification', \
		type=int, nargs='?', const=10, required=False)
	args = parser.parse_args()

	
	# check arguments are correct
	check_args(args)


	# load the data
	print('loading genetic data...')
	gtex_gct_flt = np.load(args.dataset)
	total_gene_list = np.load(args.gene_list)
	data = load_data(args.sample_json, gtex_gct_flt)


	# ensure the dataset and gene list match dimensions
	if gtex_gct_flt.shape[0] != total_gene_list.shape[0]:
		print('dataset does not match gene list.')
		sys.exit(1)


	# read subset file if provided
	if args.subset_list:
		subsets = read_subset_file(args.subset_list)
		subset_classification(data, args.gene_list, subsets, args.out_file)


	# if random is selectioned, run random 
	if args.random_test:
		random_classification(data, args.gene_list, args.num_randoms, args.rand_iters, args.out_file)



