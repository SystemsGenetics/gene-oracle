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

	# add on an additional num_per_k random samples from the data
	samples = random.sample(unused_idxs, num_per_k)
	
	for s in samples:
		final_combos.append(ast.literal_eval(sort_c_info[s][0]))

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
	return [total_gene_list[i] for i in gene_indexes]


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
	parser.add_argument('--set', help='subset to be used', type=str, required=True)
	parser.add_argument('--num_genes', help='number of genes', type=int, required=True)
	args = parser.parse_args()

	# start halo spinner
	spinner = Halo(text='Loading', spinner='dots')

	print('loading genetic data...')
	gtex_gct_flt = np.load('../datasets/gtex_gct_data_float_v7.npy')
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
			print('performing set selection via KMeans...')
			# for combos from files
			f = '../logs/' + str(args.set) + '/' + str(args.set) + '_' + str(i - 1) + '_gene_accuracy.txt'
			gene_dict = generate_new_subsets_w_clustering(f, data, total_gene_list, genes, max_experiments=80)
			# create files to write to, specify neural net architecture
			files = [str(args.set) + '_' + str(i) + '_gene_accuracy.txt']
		else:
			# for all possible combos
			gene_dict = create_raw_combos(genes, i)
			# create files to write to
			files = [str(args.set) + '_' + str(i) + '_gene_accuracy.txt']
		
		# define hidden layer sizes
		h1 = [1024]
		h2 = [1024]
		h3 = [1024]

		# open log file to write to
		fp = open('../logs/' + str(args.set) + '/' + files[0], 'w')

		for key in gene_dict:
			# retrieve the new combination of genes and create a new dataset containing the specified features
			start = time.clock()
			combo = list(key)

			gtex = GTEx(data, total_gene_list, combo)
			
			# run the neural network architecture to retrieve an accuracy based on the new dataset
			mlp = MLP(n_input=i, n_classes=53, batch_size=256, lr=0.001, epochs=75, n_h1=h1[0], n_h2=h2[0], n_h3=h3[0])
			acc = mlp.run(gtex)

			print(str(combo) + '\t' + str(acc))
			
			fp.write('{0}\t{1}\n'.format(key, acc))

		fp.close() 
