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
import json
import time

sys.path.append(os.path.dirname(os.getcwd()))

from subset_gene_test import create_random_subset, load_data, check_args, read_subset_file
from models.nn_gtex import MLP
from GTEx import GTEx


# perform random classification based on the given parameters
def random_classification(data, total_gene_list, config, num_genes, iters, out_file):
	f = open(out_file, 'w')
	f.write('Num\tAverage\tStd Dev\n')

	for num in num_genes:
		print('classifying ' + str(num) + ' random genes ' + str(iters) + ' times')
		accs = []
		# set up the neural network
		mlp = MLP(n_input=num, n_classes=len(data), batch_size=config['mlp']['batch_size'], \
			lr=config['mlp']['lr'], epochs=config['mlp']['epochs'], \
			act_funcs=config['mlp']['act_funcs'], n_layers=config['mlp']['n_h_layers'], \
			h_units=config['mlp']['n_h_units'], verbose=config['mlp']['verbose'], \
			load=config['mlp']['load'], dropout=config['mlp']['dropout'], \
			disp_step=config['mlp']['display_step'], confusion=config['mlp']['confusion'])
		
		for i in xrange(iters):
			# generate random set of genes from the total gene list
			r_genes = create_random_subset(num, total_gene_list)

			# set up the gtex class to partition data
			gtex = GTEx(data, total_gene_list, r_genes)

			# run the neural net
			accs.append(mlp.run(gtex))

		# print the results to a file
		accs_np = np.asarray(accs)
		mean = np.mean(accs_np)
		std = np.std(accs_np)
		mx = np.max(accs_np)
		mn = np.min(accs_np)
		f.write(str(num) + '\t' + str(mean) + '\t' + str(std) + '\t' + str(mx) + '\t' + str(mn) + '\n')

	f.close()


# perform classificaiton on each of the subsets provided in the subset_list argument
def subset_classification(data, total_gene_list, config, subsets, out_file, kfold_val=1):
	f = open(out_file, 'w')
	f.write('Num\tAverage\tStd Dev\tMax\tMin\n')

	for s in subsets:
		print('classifying ' + str(s))
		accs = []
		# set up the neural network
		
		for i in xrange(kfold_val):
			# set up the gtex class to partition data
			start = time.clock()
			gtex = GTEx(data, total_gene_list, subsets[s])

			mlp = MLP(n_input=gtex.train.data.shape[1], n_classes=len(data), \
				batch_size=config['mlp']['batch_size'], \
				lr=config['mlp']['lr'], epochs=config['mlp']['epochs'], \
				act_funcs=config['mlp']['act_funcs'], n_layers=config['mlp']['n_h_layers'], \
				h_units=config['mlp']['n_h_units'], verbose=config['mlp']['verbose'], \
				load=config['mlp']['load'], dropout=config['mlp']['dropout'], \
				disp_step=config['mlp']['display_step'], confusion=config['mlp']['confusion'])

			# run the neural net
			acc = mlp.run(gtex)
			accs.append(acc)
			stop = time.clock()
			print('acc is ' + str(acc) + ' time is ' + str(stop - start))

		# print the results to a file
		accs_np = np.asarray(accs)
		mean = np.mean(accs_np)
		std = np.std(accs_np)
		mx = np.max(accs_np)
		mn = np.min(accs_np)
		f.write(str(s) + '\t' + str(mean) + '\t' + str(std) + '\t' + str(mx) + '\t' + str(mn) + '\n')

	f.close()


# perform classification on every gene
def full_classification(data, total_gene_list, config, out_file, kfold_val=1):
	f = open(out_file, 'w')
	f.write('Num\tAverage\tStd Dev\n')	
	accs = []
	
	for i in xrange(kfold_val):
		# set up the gtex class to partition data
		gtex = GTEx(data, total_gene_list)

		mlp = MLP(n_input=gtex.train.data.shape[1], n_classes=len(data), \
			batch_size=config['mlp']['batch_size'], \
			lr=config['mlp']['lr'], epochs=config['mlp']['epochs'], \
			act_funcs=config['mlp']['act_funcs'], n_layers=config['mlp']['n_h_layers'], \
			h_units=config['mlp']['n_h_units'], verbose=config['mlp']['verbose'], \
			load=config['mlp']['load'], dropout=config['mlp']['dropout'], \
			disp_step=config['mlp']['display_step'], confusion=config['mlp']['confusion'])

		# run the neural net
		acc = mlp.run(gtex)
		accs.append(acc)

	# print the results to a file
	accs_np = np.asarray(accs)
	mean = np.mean(accs_np)
	std = np.std(accs_np)
	s = 'EVERY FEATURE'
	print(str(s) + '\t' + str(mean))
	f.write(str(s) + '\t' + str(mean) + '\t' + str(std) + '\n')

	f.close()



if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Run classification on specified dataset, \
		subset of genes, or a random set')
	parser.add_argument('--dataset', help='dataset to be used', type=str, required=True)
	parser.add_argument('--gene_list', help='list of genes in dataset (same order as dataset)', \
		type=str, required=True)
	parser.add_argument('--sample_json', help='json file containing number of samples per class', \
		type=str, required=True)
	parser.add_argument('--config', help='json file containing network specifications', type=str, \
		required=True)
	parser.add_argument('--out_file', help='output file to send results to', type=str, required=True)
	parser.add_argument('--subset_list', help='gmt/gct file containing subsets', type=str, required=False)
	parser.add_argument('--set', help='specific subset to run', type=str, required=False)
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

	config = json.load(open(args.config))

	# read subset file if provided
	if args.subset_list:
		subsets = read_subset_file(args.subset_list)
		
		print('checking for valid genes...')
		for s in subsets:
			genes = []
			for g in subsets[s]:
				if g in total_gene_list:
					genes.append(g)
			subsets[s] = genes
					#print('missing gene ' + str(g))
		print('done check')

		if args.set:
			sub = {}
			sub[args.set.upper()] = subsets[args.set.upper()]
			subsets = sub

		subset_classification(data, total_gene_list, config, subsets, args.out_file, kfold_val=10)


	# if random is selectioned, run random 
	if args.random_test:
		random_classification(data, total_gene_list, config, args.num_random_genes, args.rand_iters, args.out_file)


	# if not subset test and random test, run classifier on all 56k genes
	if not args.random_test and not args.subset_list:
		full_classification(data, total_gene_list, config, args.out_file)



