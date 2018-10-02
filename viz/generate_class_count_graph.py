import numpy as np 
import sys
import os
import json
from utils.utils import create_random_subset, load_data, \
					check_args, read_subset_file, \
					create_random_subset_from_interactions, \
					create_random_subset_from_NON_interactions
from models.mlp import MLP
from utils.dataset import DataContainer as DC

import matplotlib.pyplot as plt

# gtex_gct_flt = np.load('./data/float_data/gtex_gct_data_float_v7.npy')
# total_gene_list = np.load('./data/gene_lists/gtex_gene_list_v7.npy')
# data = load_data('./data/class_counts/gtex_tissue_count_v7.json', gtex_gct_flt)

gtex_gct_flt = np.load('./data/float_data/panTCGA_float_data_v2_official_symbols.npy')
total_gene_list = np.load('./data/gene_lists/gene_list_panTCGA_official_symbols.npy')
data = load_data('./data/class_counts/panTCGA_class_counts_v2.json', gtex_gct_flt)

counts = {}

for k in data:
	counts[k] = data[k].shape[1]

i = np.arange(len(counts.keys()))

colour = ["#348ABD", "#A60628"]

vals = [counts[k] for k in sorted(counts.keys())]

plt.bar(i, vals, color=colour[1], label="Sample Count", alpha=0.5, edgecolor=colour[1], lw="3")

plt.xticks(i, sorted(counts.keys()), rotation=70, ha="right")

plt.title('panTCGA Tissue Count', weight='bold')

plt.ylabel("Sample Count", weight='bold')

plt.tight_layout()

plt.show()

