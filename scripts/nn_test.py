#
# Quick test program to see if network runs
#

import sys
import os
import json
import numpy as np

sys.path.append(os.path.dirname(os.getcwd()))

from models.nn_gtex import MLP
from GTEx import GTEx

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

print('loading genetic data...')
gtex_gct_flt = np.load('../datasets/gtex_gct_data_float.npy')
total_gene_list = np.load('../datasets/gtex_complete_gene_list_str.npy')
print('done')

data = load_data("../data_scripts/numsamples.json", gtex_gct_flt)

sub = np.load('../datasets/hallmark_numpys/HALLMARK_HEDGEHOG_SIGNALING.npy')
genes = sub[:,1].tolist()

gtex = GTEx(data, total_gene_list, genes, 70, 30)

