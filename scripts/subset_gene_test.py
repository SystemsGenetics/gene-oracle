#/usr/bin/python

'''
	This file takes divides a hallmark subset into user-specified subsets in order to
	determine the best combination of genes for classification purposes.

	Prototpes:
		- run_kmeans(data,total_gene_list,combos)
		- generate_new_subsets_wo_clustering(file, data, total_gene_list, genes, max_experiments=50, rand_exps_perct=0.5)
	Todo:
		- generate_new_subsets_w_clustering
'''

import numpy as np
import os
import subprocess
import json
import shutil
import sys, argparse
import time
from halo import Halo
from sklearn.cluster import KMeans
from math import log
import ast
import random
import operator

sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.getcwd())

from models.mlp import MLP
from utils.dataset import DataContainer as DC
from utils.utils import get_combos_and_accs, load_data, check_args, \
						create_random_subset, create_raw_combos, read_subset_file

# USAGE:
# 		- perform kmeans clustering based on given paramters
# PARAMS:
#	data: dataset to be used
#	total_gene_list: list of genes in dataset (same order as dataset)
#	combos: number of combos
# NOTES:
#	- this function is not being used at the moment
def run_kmeans(data,total_gene_list,combos):
	#create data matrix of old combinations
	gene_set_data = convert_sets_to_vecs(data, total_gene_list, combos, len(combos[0]))

	inertias = []
	models = []

	# run k means k times
	print("Running Kmeans")
	for i in xrange(1,5):
		print(str(i))
		kmeans = KMeans(n_clusters=i, n_jobs=1, n_init=30, precompute_distances=False, copy_x=False)
		kmeans.fit(gene_set_data)

		models.append(kmeans)
		inertias.append(kmeans.inertia_)

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

		#extract the top num_per_k for each cluster, add to dictionary that contains tuple of classification
		#rate and cluster label

# USAGE:
# 		- perform kmeans clustering based on given paramters
# PARAMS:
#	file:file string that is an accuracy file with a list of genes
# 		 separated by a tab, followed by the accuracy for that list
#	data:
#	total_gene_list:
#	genes:
#	max_experiments:
#	rand_exps_perct
# RETURNS:
#	- returns a dictionary of new combinations with one extra gene appended that was not previously in the list
#
# NOTES:
#	-If kmeans is being run: It chooses subsets by performing KMeans
# 		clustering, choosing top performing subsets from each cluster, then also adding in some random subsets

def generate_new_subsets_wo_clustering(file, data, total_gene_list, genes, max_experiments=50, rand_exps_perct=0.5):

	# get combos and previous accuracies of the last run
	combos, prev_accs = get_combos_and_accs(file)

	combo_info = {}
	for i in xrange(len(combos)):
		combo_info[str(combos[i])] = (prev_accs[i])#, final_model.labels_[i])

	# sort the combo info descending by accuracy
	sort_c_info = sorted(combo_info.items(), key=operator.itemgetter(1), reverse=True)

	# retrieve the top num_per_k values from each cluster
	final_combos = []
	unused_idxs = []
	cnt = 0
	nxt_items = 0

	for item in sort_c_info:
		if nxt_items < max_experiments - (max_experiments * rand_exps_perct):
			final_combos.append(ast.literal_eval(item[0]))
			nxt_items = nxt_items + 1
		else:
			unused_idxs.append(cnt)
		cnt = cnt + 1

	# fill final combos with random samples from the data
	# if len(unused_idxs) > max_experiments - len(final_combos):
	# 	samples = random.sample(unused_idxs, max_experiments - len(final_combos))
	# else:
	# 	samples = unused_idxs

	if len(unused_idxs) > max_experiments * rand_exps_perct:
		samples = random.sample(unused_idxs, int(max_experiments * rand_exps_perct))
	else:
		samples = unused_idxs

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

	# get only unique combinations
	for i in next_set_size_combos:
		i.sort()
	unique = [list(x) for x in set(tuple(x) for x in next_set_size_combos)]

	ret_combos = []
	for f in unique:
		ret_combos.append(tuple(f))

	return dict.fromkeys(ret_combos)


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
	parser.add_argument('--num_genes', help='number of genes', type=int, required=False)
	parser.add_argument('--log_dir', help='directory where logs are stored', type=str, required=True)
	args = parser.parse_args()

	# check arguments are correct
	check_args(args)

	# start halo spinner
	spinner = Halo(text='Loading', spinner='dots')

	print('loading genetic data...')
	gtex_gct_flt = np.load(args.dataset)
	total_gene_list = np.load(args.gene_list, encoding='ASCII')
	print('done')

	data = load_data(args.sample_json, gtex_gct_flt)

	# load the data
	if "random" in args.set:
		genes = create_random_subset(args.num_genes, total_gene_list)

	else:
		if args.subset_list:
			subsets = read_subset_file(args.subset_list)
			for s in subsets:
				genes = []
				for g in subsets[s]:
					if g in total_gene_list:
						genes.append(g)
				subsets[s] = genes

			try:
				genes = subsets[args.set.upper()]
			except:
				print('Set not found in subset file, try again')
				sys.exit(1)
		else:
			print('must include subset file if not performing random test. exiting.')
			sys.exit(1)

	config = json.load(open(args.config))

	if not os.path.exists(args.log_dir):
		os.makedirs(args.log_dir)

	with open(args.log_dir + '/gene_list.txt', 'w') as f:
		for i in genes:
			f.write(str(i) + '\n')
		f.close()

	print('beginning search for optimal combinations...')
	for i in xrange(1, len(genes) + 1):
		print('--------ITERATION ' + str(i) + '--------')

		# read in the previous accuracy file
		if i > 3 and i != args.num_genes:
			# for combos from files
			f = args.log_dir + '/' + str(args.set) + '_' + str(i - 1) + '_gene_accuracy.txt'
			gene_dict = generate_new_subsets_wo_clustering(f, data, total_gene_list, genes, \
				max_experiments=60, rand_exps_perct=0.5)
		else:
			#for all possible combos
			gene_dict = create_raw_combos(genes, i)

		files = [str(args.set) + '_' + str(i) + '_gene_accuracy.txt']

		# open log file to write to
		fp = open(args.log_dir + '/' + files[0], 'w')

		print("Running MLP")
		mlp = MLP(n_input=i, n_classes=len(data), \
				batch_size=config['mlp']['batch_size'], \
				lr=config['mlp']['lr'], epochs=config['mlp']['epochs'], \
				act_funcs=config['mlp']['act_funcs'], n_layers=config['mlp']['n_h_layers'], \
				h_units=config['mlp']['n_h_units'], verbose=config['mlp']['verbose'], \
				load=config['mlp']['load'], dropout=config['mlp']['dropout'], \
				disp_step=config['mlp']['display_step'], confusion=config['mlp']['confusion'])

		for key in gene_dict:
			# retrieve the new combination of genes and create a new dataset containing the specified features
			#start = time.clock()
			combo = list(key)

			dataset = DC(data, total_gene_list, combo)
			#stop = time.clock()
			#print('data load: ' + str(stop - start))
			start = time.clock()
			# run the neural network architecture to retrieve an accuracy based on the new dataset
			acc = mlp.run(dataset)
			stop = time.clock()

			#print('time nn: ' + str(stop - start))
			print(str(combo) + '\t' + str(acc) + '\t' + str(stop - start))

			fp.write('{0}\t{1}\n'.format(str(key), acc))

		fp.close()
