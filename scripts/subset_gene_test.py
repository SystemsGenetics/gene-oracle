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
from halo import Halo
from sklearn.cluster import KMeans
from math import log
import ast
import random
import operator 
import re
import json


sys.path.append(os.path.dirname(os.getcwd()))

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

# special helper function to sanitize the string containing the genes from
# an accuracy file
def sanitize(gene_str):
	gene_list = gene_str.strip('()')
	gene_list = gene_list.replace('\'', '')
	gene_list = gene_list.replace(' ', '')
	gene_list = gene_list.split(',')
	return gene_list

# comment this yo
# convert sets 
def convert_sets_to_vecs(data, total_gene_list, combo_list, set_size):
	feature_list = []
	for combo in combo_list:
		dataset = GTEx(data, total_gene_list, combo, train_split=30, test_split=70)

		concat_genes = dataset.train.data[:,0]

		for i in xrange(1, set_size):
			concat_genes = np.append(concat_genes, dataset.train.data[:,i])

		feature_list.append(concat_genes)

	# convert to numpy format
	x_data = np.array(feature_list)

	return x_data

# generate_new_subsets_w_clustering takes in a file string that is an accuracy file with a list of genes
# separated by a tab, followed by the accuracy for that list. it returns a dictionary of new combinations 
# with one extra gene appended that was not previously in the list. It chooses subsets by performing KMeans
# clustering, choosing top performing subsets from each cluster, then also adding in some random subsets
def generate_new_subsets_w_clustering(file, data, total_gene_list, genes, max_experiments=80):
	# collect previous files combinations/accuracyies
	prev_combos = []
	prev_run = np.loadtxt(file, delimiter='\t', dtype=np.str)
	
	# gather previous combinations
	combos = []
	prev_combos = prev_run[:,0]
	for pc in prev_combos:
		combos.append(sanitize(pc))

	# gather previous accuracies
	prev_accs = prev_run[:,1]

	# create data matrix of old combinations
	gene_set_data = convert_sets_to_vecs(data, total_gene_list, combos, len(combos[0]))

	inertias = []
	BIC_list = []
	models = []

	# run k means k times
	print("Running Kmeans")
	for i in xrange(1,11):
		kmeans = KMeans(n_clusters=i, n_jobs=8, n_init=30)
		kmeans.fit(gene_set_data)

		models.append(kmeans)
		inertias.append(kmeans.inertia_)

		# calculate BIC and append to list
		BIC = log(kmeans.inertia_) - (log(len(combos) * i))
		BIC_list.append(BIC)

	print('done kmeans')

	# approximate second derivatives to determine where the 'elbow' in the curve is
	second_dervs = []
	for i in xrange(1, len(inertias) - 1):
		xpp = inertias[i + 1] + inertias[i - 1] - 2 * inertias[i]
		second_dervs.append(xpp)

	# add one... excluded first and last k from calculations TODO: may need to fix this
	final_k = second_dervs.index(max(second_dervs)) + 1
	final_model = models[final_k]

	print('final k for ' + str(len(combos[0])) + ' combos is ' + str(final_k))

	# find the top num sets from each cluster and additionally return num random sets
	# num = max_experiments / (k + 1) send off num sets from each k clusters + num random sets
	num_per_k = max_experiments / (final_k + 2)

	# extract the top num_per_k for each cluster, add to dictionary that contains tuple of classification 
	# rate and cluster label
	combo_info = {}
	for i in xrange(len(combos)):
		combo_info[str(combos[i])] = (prev_accs[i], final_model.labels_[i])

	# sort the combo info descending by accuracy
	sort_c_info = sorted(combo_info.items(), key=operator.itemgetter(1), reverse=True)

	# retrieve the top num_per_k values from each cluster
	final_combos = []
	unused_idxs = []
	cnt = 0
	counts = np.zeros(final_k + 1)
	for item in sort_c_info:
		if counts[item[1][1]] < num_per_k: # see if the num_per_k is under threshold 
			final_combos.append(ast.literal_eval(item[0])) # append the list of genes to final_combos
			counts[item[1][1]] = counts[item[1][1]] + 1 # increment the global counts of each cluster
		else:
			unused_idxs.append(cnt)		
		
		cnt = cnt + 1 # increment index count

	# fill final combos with random samples from the data
	samples = random.sample(unused_idxs, max_experiments - len(final_combos))
	
	for s in samples:
		final_combos.append(ast.literal_eval(sort_c_info[s][0]))

	print('num sets moving on is ' + str(len(final_combos)))

	# append missing genes to each of the combinations of 3
	next_set_size_combos = []
	for c in final_combos:
		for g in genes:
			if g not in c:
				temp_list = c[:]
				temp_list.append(g)
				next_set_size_combos.append(temp_list)

	ret_combos = []
	for f in next_set_size_combos:
		ret_combos.append(tuple(f))

	return dict.fromkeys(ret_combos)

# create every possible combination
def create_raw_combos(genes, i):
	combos = []
	for c in itertools.combinations(genes, i):
		combos.append(c)

	return dict.fromkeys(combos)

# get random gene indexes between 0-56238
def create_random_subset(num_genes, total_gene_list):		
	#Generate Gene Indexes for Random Sample
	gene_indexes = np.random.randint(0, len(total_gene_list), num_genes)
	return [total_gene_list[i][0] for i in gene_indexes]


def load_data(num_samples_json, gtex_gct_flt):
	sample_count_dict = {}
	with open(num_samples_json) as f:
		sample_count_dict = json.load(f)

	idx = 0
	data = {}

	for k in sorted(sample_count_dict.keys()):
		data[k] = gtex_gct_flt[:,idx:(idx + int(sample_count_dict[k]))]
		idx = idx + int(sample_count_dict[k])

	return data


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Run tests on specified subsets of a hallmark or random set')
	parser.add_argument('--dataset', help='dataset to be used', type=str, required=True)
	parser.add_argument('--gene_list', help='list of genes in dataset (same order as dataset)', \
		type=str, required=True)
	parser.add_argument('--sample_json', help='json file containing number of samples per class', \
		type=str, required=True)
	parser.add_argument('--config', help='json file containing network specifications', type=str, \
		required=True)
	parser.add_argument('--subset_list', help='gmt/gct/txt file containing subsets', type=str, required=False)
	parser.add_argument('--set', help='subset to be used', type=str, required=True)
	parser.add_argument('--num_genes', help='number of genes', type=int, required=True)
	parser.add_argument('--log_dir', help='directory where logs are stored', type=str, required=True)
	args = parser.parse_args()

	# check arguments are correct
	check_args(args)

	# start halo spinner
	spinner = Halo(text='Loading', spinner='dots')

	print('loading genetic data...')
	gtex_gct_flt = np.load(args.dataset)
	total_gene_list = np.load(args.gene_list)
	print('done')

	data = load_data(args.sample_json, gtex_gct_flt)

	# load the hedgehog data
	if "random" in args.set:
		genes = create_random_subset(args.num_genes, total_gene_list)
	else:
		if args.subset_list:
			subsets = read_subset_file(args.subset_list)
			try:
				genes = subsets[args.set.upper()]
				for g in genes:
					if g not in total_gene_list:
						genes.remove(g)
						print('missing gene ' + str(g))
			except:
				print('set not found in subset file, try again')
				sys.exit(1)

		else:
			print('must include subset file if not performing random test. exiting.')
			sys.exit(1)
		
	config = json.load(open(args.config))

	if not os.path.exists(args.log_dir):
		os.makedirs(args.log_dir)

	print('beginning search for optimal combinations...')
	for i in xrange(6, len(genes)):
		print('--------ITERATION ' + str(i) + '--------')

		# read in the previous accuracy file
		if i > 3:
			print('performing set selection via KMeans...')
			# for combos from files
			f = args.log_dir + str(args.set) + '_' + str(i - 1) + '_gene_accuracy.txt'
			gene_dict = generate_new_subsets_w_clustering(f, data, total_gene_list, genes, max_experiments=80)
		else:
			# for all possible combos
			gene_dict = create_raw_combos(genes, i)
		
		files = [str(args.set) + '_' + str(i) + '_gene_accuracy.txt']

		# open log file to write to
		fp = open(args.log_dir + '/' + files[0], 'w')

		mlp = MLP(n_input=i, n_classes=len(data), \
				batch_size=config['mlp']['batch_size'], \
				lr=config['mlp']['lr'], epochs=config['mlp']['epochs'], \
				act_funcs=config['mlp']['act_funcs'], n_layers=config['mlp']['n_h_layers'], \
				h_units=config['mlp']['n_h_units'], verbose=config['mlp']['verbose'], \
				load=config['mlp']['load'], dropout=config['mlp']['dropout'], \
				disp_step=config['mlp']['display_step'], confusion=config['mlp']['confusion'])

		for key in gene_dict:
			# retrieve the new combination of genes and create a new dataset containing the specified features
			start = time.clock()
			combo = list(key)

			gtex = GTEx(data, total_gene_list, combo)
			
			# run the neural network architecture to retrieve an accuracy based on the new dataset
			acc = mlp.run(gtex)

			print(str(combo) + '\t' + str(acc))
			
			fp.write('{0}\t{1}\n'.format(key, acc))

		fp.close() 
