import json
import re
import itertools
import os
import numpy as np
import sys


# USAGE:
#	-check the arguments are correct for the program
# PARAMS:
#	args: list of program arguments
def check_args(args):
	# check dataset is of correct type
	if os.path.exists(args.dataset):
		split = args.dataset.split('.')
		if split[-1] != 'npy':
			print('Dataset file must be a numpy file.')
			sys.exit(1)
	else:
		print('File ' + str(args.dataset) + ' does not exist!')
		sys.exit(1)

	# check gene list is of correct type
	if os.path.exists(args.gene_list):
		split = args.gene_list.split('.')
		if split[-1] != 'npy':
			print('Gene list file must be a numpy file.')
			sys.exit(1)
	else:
		print('File ' + str(args.gene_list) + ' does not exist!')
		sys.exit(1)

	# check gene list is of correct type
	if os.path.exists(args.sample_json):
		split = args.sample_json.split('.')
		if split[-1] != 'json':
			print('sample file must be a json file.')
			sys.exit(1)
	else:
		print('File ' + str(args.sample_json) + ' does not exist!')
		sys.exit(1)

# USAGE:
# 	- create every possible combination
# PARAMS:
#	genes:
#	i:
# Todo: rename i to something better....@Colin
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

# USAGE:
# 	- read a csv or txt file that contains a name of a subset followed by a list of genes
# PARAMS:
#	file: file to read
def read_subset_file(file):
	with open(file, 'r') as f:
		content = f.readlines()

	# eliminate new line characters
	content = [x.strip() for x in content]

	# split on tabs or commas to create a sublist of set names and genes
	content = [re.split('\t|,| ', x) for x in content]

	# create a dictionary with keys subset names and values list of genes
	subsets = {}
	for c in content:
		subsets[c[0]] = c[1:]

	return subsets



# USAGE:
#  -special helper function to sanitize the string containing the genes from
# an accuracy file
# PARAMS:
#	gene_str: gene string to sanitize
def sanitize(gene_str):
	gene_list = gene_str.strip('()')
	gene_list = gene_list.replace('\'', '')
	gene_list = gene_list.replace(' ', '')
	gene_list = gene_list.split(',')
	return gene_list


# USAGE:
# - convert sets
# PARAMS:
#	data:
#	total_gene_list:
# 	combo_list:
#	set_size:
def convert_sets_to_vecs(data, total_gene_list, combo_list, set_size):
	feature_list = []
	for combo in combo_list:
		dataset = GTEx(data, total_gene_list, combo, train_split=100, test_split=0)

		concat_genes = dataset.train.data[:,0]

		for i in xrange(1, set_size):
			concat_genes = np.append(concat_genes, dataset.train.data[:,i])

		feature_list.append(concat_genes)

	# convert to numpy format
	x_data = np.array(feature_list)

	return x_data


# USAGE:
# - return the combinations and accuracies that are listed in a log file
# PARAMS:
#	file:
def get_combos_and_accs(file):
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

	return combos, prev_accs


# USAGE:
# - return a list of subsets only containing genes with interactions within the subset
# PARAMS:
# subsets - dictionary of subsets similar to whats returned by read_subset_file
# interaction_file - file containing gene interactions

# NOTE: the interaction file should not contain any metadata in the top (i.e. a readme)...
# it should just contain the column headers followed by the interaction data.
# column format is expected to follow the biogrid interaction file 
def convert_subset_to_interactions(subsets, interaction_file):
	# read biogrid file into a dataframe, then extract official symbols
	biogrid = pd.read_csv(interaction_file, sep='\t')
	interactions = biogrid[['OFFICIAL_SYMBOL_A', 'OFFICIAL_SYMBOL_B']]
	interactions = interactions.values
	interactions = interactions.astype(np.str)

	interacted_subsets = {}
	for s in subsets:
		interacted_subsets[s] = []
		subset_interacting_genes = []
		
		for g in subsets[s]:
			total_interacting_genes = []
			locs = np.where(interactions==g)
			
			# go through each interacting gene and add it
			for i in range(locs[0].shape[0]):
				# doing 1 - locs[1][i] gives you the interacting gene
				cand_gene = interactions[locs[0][i], 1 - locs[1][i]]

				if cand_gene not in total_interacting_genes:
					total_interacting_genes.append(cand_gene)

			for gene in total_interacting_genes:
				if gene in subsets[s] and gene not in interacted_subsets[s]:
					interacted_subsets[s].append(gene)

	return interacted_subsets








