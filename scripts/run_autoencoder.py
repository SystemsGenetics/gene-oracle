# 
# temp script for running autoencoder
#

import numpy as np
import sys, os, argparse

sys.path.append(os.path.dirname(os.getcwd()))

from subset_gene_test import convert_sets_to_vecs, load_data, get_combos_and_accs
from models.autoencoder import autoencoder





if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Run tests on specified subsets of a hallmark or random set')
	parser.add_argument('--dataset', help='dataset to be used', type=str, required=True)
	parser.add_argument('--gene_list', help='list of genes in dataset (same order as dataset)', \
		type=str, required=True)
	parser.add_argument('--sample_json', help='json file containing number of samples per class', \
		type=str, required=True)
	parser.add_argument('--file', help='prev file to read from', type=str, required=True)
	args = parser.parse_args()

	# load genetic data and the list of genes associated
	print('loading genetic data...')
	gtex_gct_flt = np.load(args.dataset)
	total_gene_list = np.load(args.gene_list)
	print('done')

	# load data into dictionary for easy access
	data = load_data(args.sample_json, gtex_gct_flt)

	# get combos and accs of last run
	combos, _ = get_combos_and_accs(args.file)

	# gather data into 
	# print('loading data for autoencoder')
	# x_data = convert_sets_to_vecs(data, total_gene_list, combos, len(combos[0]))
	# print(x_data.shape)

	x_data = np.load('./temp_autoencoder_data.npy')	

	print('training autoencoder')
	ace = autoencoder(n_input=x_data.shape[1])
	ace.run(x_data)
