#/usr/bin/python

'''
	This script can be used to run a specified dataset, a specified subset of genes,
	or a specified number of random genes for classification.

	It is required to have a numpy array containing column-wise samples and rows of
	genes. Additionally, it is required to have a numpy vector of genes that are contained
	in the dataset (in the same exact order).


	Protypes:
	- random_classification(data, total_gene_list, config, num_genes, iters, out_file, kfold_val)
	- subset_classification(data, total_gene_list, config, subsets, out_file, kfold_val=1)
	- full_classification(data, total_gene_list, config, out_file, kfold_val=1)

	Todo:
		-modularize functions in main, main should be small, logic needs to be in functions

'''

import numpy as np
import sys, argparse
import os
import json
import time

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())

from utils.utils import create_random_subset, load_data, check_args, read_subset_file
from models.mlp import MLP
from utils.dataset import DataContainer as DC


# USAGE:
# 		- perform random classification based on the given parameters
# PARAMS:
#	data: dataset to be used
#	total_gene_list: list of genes in dataset (same order as dataset)
#	config: json file containing network specifications
#	num_genes: Number of random genes to assess
#	iters: Number of iterations to perform for random classification
#	out_file: the file to write to
#	kfold_val: Number of folds for K fold cross validation
def random_classification(data, total_gene_list, config, num_genes, iters, out_file, kfold_val):
	if out_file:
		f = open(out_file, 'w')
		f.write('Num\tAverage\tStd Dev\n')

	for num in num_genes:
		print('classifying ' + str(num) + ' random genes ' + str(iters) + ' times')
		# set up the neural network
		mlp = MLP(n_input=num, n_classes=len(data), batch_size=config['mlp']['batch_size'], \
			lr=config['mlp']['lr'], epochs=config['mlp']['epochs'], \
			act_funcs=config['mlp']['act_funcs'], n_layers=config['mlp']['n_h_layers'], \
			h_units=config['mlp']['n_h_units'], verbose=config['mlp']['verbose'], \
			load=config['mlp']['load'], dropout=config['mlp']['dropout'], \
			disp_step=config['mlp']['display_step'])#, confusion=config['mlp']['confusion'], roc=config['mlp']['roc'], pr=config['mlp']['pr'])

		for i in xrange(iters):
			# generate random set of genes from the total gene list
			r_genes = create_random_subset(num, total_gene_list)
			accs = []
			for _ in xrange(kfold_val):
				# set up the DataContainer class to partition data
				dataset = DC(data, total_gene_list, r_genes)

				# run the neural net
				accs.append(mlp.run(dataset))

			# print the results to a file
			accs_np = np.asarray(accs)
			mean = np.mean(accs_np)
			std = np.std(accs_np)
			mx = np.max(accs_np)
			mn = np.min(accs_np)
			print(str(num) + '\t' + str(mean) + '\t' + str(std) + '\t' + str(mx) + '\t' + str(mn))
			if out_file:
				f.write(str(num) + '\t' + str(mean) + '\t' + str(std) + '\t' + str(mx) + '\t' + str(mn) + '\n')

	if out_file:
		f.close()

# USAGE:
# 		- perform classification on each of the subsets provided in the subset_list argument
# PARAMS:
#	data: dataset to be used
#	total_gene_list: list of genes in dataset (same order as dataset)
#	config: json file containing network specifications
#	subsets: gmt/gct file containing subsets
#	out_file: the file to write to
#	kfold_val: Number of folds for K fold cross validation, set at 1
def subset_classification(data, total_gene_list, config, subsets, out_file, kfold_val=1):
	if out_file:
		f = open(out_file, 'w')
		f.write('Num\tAverage\tStd_Dev\tMax\tMin\n')

	for s in subsets:
		accs = []
		# set up the neural network

		for i in xrange(kfold_val):
			# set up the gtex class to partition data
			dataset = DC(data, total_gene_list, subsets[s])

			mlp = MLP(n_input=dataset.train.data.shape[1], n_classes=len(data), \
				batch_size=config['mlp']['batch_size'], \
				lr=config['mlp']['lr'], epochs=config['mlp']['epochs'], \
				act_funcs=config['mlp']['act_funcs'], n_layers=config['mlp']['n_h_layers'], \
				h_units=config['mlp']['n_h_units'], verbose=config['mlp']['verbose'], \
				load=config['mlp']['load'], dropout=config['mlp']['dropout'])#, \
				#disp_step=config['mlp']['display_step'], confusion=config['mlp']['confusion'], roc=config['mlp']['roc'],pr=config['mlp']['pr'])

			# run the neural net
			acc = mlp.run(dataset)
			accs.append(acc)
			#print('acc is ' + str(acc) + ' time is ' + str(stop - start))

		# print the results to a file
		accs_np = np.asarray(accs)
		mean = np.mean(accs_np)
		std = np.std(accs_np)
		mx = np.max(accs_np)
		mn = np.min(accs_np)
		print(str(s) + '\t' + str(mean) + '\t' + str(std) + '\t' + str(mx) + '\t' + str(mn))
		if out_file:
			f.write(str(s) + '\t' + str(mean) + '\t' + str(std) + '\t' + str(mx) + '\t' + str(mn) + '\n')

	if out_file:
		f.close()

# USAGE:
# 		- perform classification on every gene
# PARAMS:
#	data: dataset to be used
#	total_gene_list: list of genes in dataset (same order as dataset)
#	config: json file containing network specifications
#	out_file: the file to write to
#	kfold_val: Number of folds for K fold cross validation, set at 1
def full_classification(data, total_gene_list, config, out_file, kfold_val=1):
	if out_file:
		f = open(out_file, 'w')
		f.write('Num\tAverage\tStd Dev\n')

	accs = []

	for i in xrange(kfold_val):
		# set up the data class to partition data
		dataset = DC(data, total_gene_list)

		mlp = MLP(n_input=dataset.train.data.shape[1], n_classes=len(data), \
			batch_size=config['mlp']['batch_size'], \
			lr=config['mlp']['lr'], epochs=config['mlp']['epochs'], \
			act_funcs=config['mlp']['act_funcs'], n_layers=config['mlp']['n_h_layers'], \
			h_units=config['mlp']['n_h_units'], verbose=config['mlp']['verbose'], \
			load=config['mlp']['load'], dropout=config['mlp']['dropout'], \
			disp_step=config['mlp']['display_step'], confusion=config['mlp']['confusion'], roc=config['mlp']['roc'],pr=config['mlp']['pr'])

		# run the neural net
		acc = mlp.run(dataset)
		accs.append(acc)

	# print the results to a file
	accs_np = np.asarray(accs)
	mean = np.mean(accs_np)
	std = np.std(accs_np)
	s = 'EVERY FEATURE'
	print(str(s) + '\t' + str(mean))
	if out_file:
		f.write(str(s) + '\t' + str(mean) + '\t' + str(std) + '\n')
		f.close()



if __name__ == '__main__':

	#Parse Arguments
	parser = argparse.ArgumentParser(description='Run classification on specified dataset, \
		subset of genes, or a random set')
	parser.add_argument('--dataset', help='dataset to be used', type=str, required=True)
	parser.add_argument('--gene_list', help='list of genes in dataset (same order as dataset)', \
		type=str, required=True)
	parser.add_argument('--sample_json', help='json file containing number of samples per class', \
		type=str, required=True)
	parser.add_argument('--config', help='json file containing network specifications', type=str, \
		required=True)
	parser.add_argument('--out_file', help='output file to send results to', type=str, required=False)
	parser.add_argument('--subset_list', help='gmt/gct file containing subsets', type=str, required=False)
	parser.add_argument('--set', help='specific subset to run', type=str, required=False)
	parser.add_argument('--random_test', help='Perform random test', action='store_true', required=False)
	parser.add_argument('--num_random_genes', help='Number of random genes to assess', nargs='+', \
		type=int, required=False)
	parser.add_argument('--rand_iters', help='Number of iterations to perform for random classification', \
		type=int, nargs='?', const=10, required=False)
	parser.add_argument('--k_fold', help='Number of folds for K fold cross validation', \
		type=int, nargs='?', const=10, required=False)


	args = parser.parse_args()

	# Check arguments are correct
	check_args(args)

	# load the data
	print('loading genetic data...')
	gtex_gct_flt = np.load(args.dataset)
	total_gene_list = np.load(args.gene_list)
	data = load_data(args.sample_json, gtex_gct_flt)

	# ensure the dataset and gene list match dimensions
	#assert gtex_gct_flt.shape[0] not total_gene_list.shape[0], "dataset does not match gene list."
	if gtex_gct_flt.shape[0] != total_gene_list.shape[0]:
		print('dataset does not match gene list.')
		sys.exit(1)

	config = json.load(open(args.config))

	# RUN SUBSET CLASSIFICATION
	# read subset file, if provided
	if args.subset_list and not args.random_test:
		subsets = read_subset_file(args.subset_list)

		tot_genes = []
		missing_genes = []

		print('checking for valid genes...')
		for s in subsets:
			genes = []
			for g in subsets[s]:
				if g not in tot_genes:
					tot_genes.append(g)
				if g in total_gene_list:
					genes.append(g)
				else:
					if g not in missing_genes:
						missing_genes.append(g)
			subsets[s] = genes
					#print('missing gene ' + str(g))
		print('missing ' + str(len(missing_genes)) + '/' + str(len(tot_genes)) + ' genes' + ' or ' \
			 + str(int((float(len(missing_genes)) / len(tot_genes)) * 100.0)) + '% of genes')

		if args.set:
			sub = {}
			sub[args.set.upper()] = subsets[args.set.upper()]
			subsets = sub

		subset_classification(data, total_gene_list, config, subsets, args.out_file, kfold_val=10)


	#RUN RANDOM CLASSIFICATION
	# if random is selectioned, run random
	if args.random_test:
		if args.num_random_genes:
			random_classification(data, total_gene_list, config, args.num_random_genes, args.rand_iters, args.out_file, kfold_val=10)
		elif args.subset_list:
			# get the number of genes for each subset
			num = []
			subsets = read_subset_file(args.subset_list)
			for s in subsets:
				genes = []
				for g in subsets[s]:
					if g in total_gene_list:
						genes.append(g)
				subsets[s] = genes

			for k in subsets:
				num.append(len(subsets[k]))
			num.sort()
			random_classification(data, total_gene_list, config, num, args.rand_iters, args.out_file, kfold_val=1)


	#RUN FULL_CLASSIFICATION
	# if not subset test and random test, run classifier on all 56k genes
	if not args.random_test and not args.subset_list:
		full_classification(data, total_gene_list, config, args.out_file)
