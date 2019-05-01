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

from utils.utils import create_random_subset, load_data, \
						check_args, read_subset_file, \
						create_random_subset_from_interactions, \
						create_random_subset_from_NON_interactions
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
def random_classification(data, total_gene_list, config, \
							num_genes, iters, out_file, kfold_val, verbose=True,\
							interaction_genes=None, interaction_list=None):
	if out_file:
		f = open(out_file, 'w')
		f.write('\t'.join(['Name'] + ['%d' % (i) for i in range(kfold_val * iters)]) + '\n')

	for num in num_genes:
		if verbose:
			print('classifying ' + str(num) + ' random genes ' + str(iters) + ' times')

		# set up the neural network
		mlp = MLP(n_input=num, n_classes=len(data), batch_size=config['mlp']['batch_size'], \
			lr=config['mlp']['lr'], epochs=config['mlp']['epochs'], \
			act_funcs=config['mlp']['act_funcs'], n_layers=config['mlp']['n_h_layers'], \
			h_units=config['mlp']['n_h_units'], verbose=config['mlp']['verbose'], \
			load=config['mlp']['load'], dropout=config['mlp']['dropout'], \
			disp_step=config['mlp']['display_step'])#, confusion=config['mlp']['confusion'], roc=config['mlp']['roc'], pr=config['mlp']['pr'])

		accuracies = []

		for i in range(iters):
			# generate random set of genes from the total gene list
			if interaction_genes:
				r_genes = create_random_subset_from_NON_interactions(num, \
																total_gene_list, \
																interaction_genes, \
																interaction_list)
			else:
				r_genes = create_random_subset(num, total_gene_list)


			# set up the DataContainer class to partition data
			dataset = DC(data, total_gene_list, r_genes)

			for _ in range(kfold_val):
				# run the neural net
				accuracies.append(mlp.run(dataset))

		# print the results to a file
		line = '\t'.join(['%d' % (num)] + ['%.6f' % (acc) for acc in accuracies])

		if verbose:
			print(line)
		if out_file:
			f.write(line + '\n')

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
def subset_classification(data, total_gene_list, config, subsets, out_file, kfold_val=10, verbose=True):
	if out_file:
		f = open(out_file, 'w')
		f.write('\t'.join(['Name'] + ['%d' % (i) for i in range(kfold_val)]) + '\n')

	for subset_name in subsets:
		accuracies = []

		for i in range(kfold_val):
			# set up the gtex class to partition data
			dataset = DC(data, total_gene_list, subsets[subset_name])

			# set up the neural network
			mlp = MLP(n_input=dataset.train.data.shape[1], n_classes=len(data), \
				batch_size=config['mlp']['batch_size'], \
				lr=config['mlp']['lr'], epochs=config['mlp']['epochs'], \
				act_funcs=config['mlp']['act_funcs'], n_layers=config['mlp']['n_h_layers'], \
				h_units=config['mlp']['n_h_units'], verbose=config['mlp']['verbose'], \
				load=config['mlp']['load'], dropout=config['mlp']['dropout'])#, \
				#disp_step=config['mlp']['display_step'], confusion=config['mlp']['confusion'], roc=config['mlp']['roc'],pr=config['mlp']['pr'])

			# run the neural net
			acc = mlp.run(dataset)
			accuracies.append(acc)
			#print('acc is ' + str(acc) + ' time is ' + str(stop - start))

		# print the results to a file
		line = '\t'.join([subset_name] + ['%.6f' % (acc) for acc in accuracies])

		if verbose:
			print(line)
		if out_file:
			f.write(line + '\n')

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
def full_classification(data, total_gene_list, config, out_file, kfold_val=10, verbose=True):
	if out_file:
		f = open(out_file, 'w')
		f.write('\t'.join(['Name'] + ['%d' % (i) for i in range(kfold_val)]) + '\n')

	accuracies = []

	for i in range(kfold_val):
		# set up the data class to partition data
		dataset = DC(data, total_gene_list)

		mlp = MLP(n_input=dataset.train.data.shape[1], n_classes=len(data), \
			batch_size=config['mlp']['batch_size'], \
			lr=config['mlp']['lr'], epochs=config['mlp']['epochs'], \
			act_funcs=config['mlp']['act_funcs'], n_layers=config['mlp']['n_h_layers'], \
			h_units=config['mlp']['n_h_units'], verbose=config['mlp']['verbose'], \
			load=config['mlp']['load'], dropout=config['mlp']['dropout'], )#\
			#disp_step=config['mlp']['display_step'], confusion=config['mlp']['confusion'], roc=config['mlp']['roc'],pr=config['mlp']['pr'])

		# run the neural net
		acc = mlp.run(dataset)
		accuracies.append(acc)

	# print the results to a file
	line = '\t'.join(['EVERY_FEATURE'] + ['%.6f' % (acc) for acc in accuracies])

	if verbose:
		print(line)

	if out_file:
		f.write(line + '\n')
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
	parser.add_argument('--range_random_genes', help='range of random genes to assess', nargs='+', \
		type=int, required=False)
	parser.add_argument('--rand_iters', help='Number of iterations to perform for random classification', \
		type=int, default=50, required=False)
	parser.add_argument('--k_fold', help='Number of folds for K fold cross validation', \
		type=int, default=10, required=False)
	parser.add_argument('--verbose', help='verbose if true', action='store_true', required=False)
	parser.add_argument('--interaction_genes', help='list of valid genes from interaction list', \
		type=str, required=False, default=None)
	parser.add_argument('--interaction_list', help='pairwise list of interacting genes', \
		type=str, required=False, default=None)

	args = parser.parse_args()

	# Check arguments are correct
	check_args(args)

	# load the data
	print('loading genetic data...')
	gtex_gct_flt = np.load(args.dataset)
	total_gene_list = np.load(args.gene_list)
	data = load_data(args.sample_json, gtex_gct_flt)

	# load interaction data, if passed
	if args.interaction_genes:
		interaction_genes = np.load(args.interaction_genes)

		# ensure only genes in interaction_genes are contained within the dataset
		interaction_genes = [g for g in interaction_genes if g in total_gene_list]
	else:
		interaction_genes = None

	if args.interaction_list:
		interaction_list = np.load(args.interaction_list)
		original_interaction_genes = np.load(args.interaction_genes)

		# find missing genes, then delete them from interaction list
		missing = [g for g in original_interaction_genes if g not in total_gene_list]
		for g in missing:
			locs = np.where(interaction_list==g)
			interaction_list = np.delete(interaction_list, locs[0], axis=0)

		interaction_genes = list(np.unique(interaction_list))
	else:
		interaction_list = None

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

		subset_classification(data, total_gene_list, config, subsets, args.out_file, kfold_val=args.k_fold, verbose=args.verbose)


	#RUN RANDOM CLASSIFICATION
	# if random is selectioned, run random
	if args.random_test:
		if args.range_random_genes:
			gene_range = range(args.range_random_genes[0], args.range_random_genes[1] + 1)
			random_classification(data, total_gene_list, config, gene_range, \
									args.rand_iters, args.out_file, kfold_val=args.k_fold, verbose=args.verbose)
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
			random_classification(data, total_gene_list, config, num, args.rand_iters, \
									args.out_file, kfold_val=3, interaction_genes=interaction_genes, \
									interaction_list=interaction_list)


	#RUN FULL_CLASSIFICATION
	# if not subset test and random test, run classifier on all 56k genes
	if not args.random_test and not args.subset_list:
		full_classification(data, total_gene_list, config, args.out_file)
