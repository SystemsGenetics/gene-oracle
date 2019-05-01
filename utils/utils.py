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

# get random gene indexes between 0-len total_gene_list
def create_random_subset(num_genes, total_gene_list):
	#Generate Gene Indexes for Random Sample
	gene_indexes = np.random.randint(0, len(total_gene_list), num_genes)
	return [total_gene_list[i] for i in gene_indexes]


# get random genes that are within an interaction list
# num_genes: number of random genes to select
# total_gene_list: list of genes within the dataset being used
# interaction_genes: list of valid genes from an interaction list
# interaction_list: should be two columns: interacting A and interacting B...
#					a list of interacting genes (pairwise)
# NOTE: it is assumed that interaction_genes are only genes contained inside of
# total_gene_list, perform a check before passing it to this function
def create_random_subset_from_interactions(num_genes, total_gene_list, \
											interaction_genes, interaction_list):
	num_interactions = num_genes / 2
	set_length = 0

	while set_length != num_interactions * 2:
		gene_indexes = np.random.randint(0, len(interaction_genes), num_interactions)
		rand_genes = [interaction_genes[i] for i in gene_indexes]

		final_genes = []
		for g in rand_genes:
			# get locations for g in the interaction list
			locs = np.where(interaction_list == g)
			
			# get a random interaction within the available locations
			if locs[0].shape[0] > 1:
				idx = np.random.randint(0, locs[0].shape[0], 1)[0]
			else:
				idx = 0

			final_genes.append(interaction_list[locs[0][idx]][0])
			final_genes.append(interaction_list[locs[0][idx]][1])

		# if odd amount of genes, add one extra gene
		if num_genes % 2:
			idx = np.random.randint(0, len(interaction_genes), 1)
			final_genes.append(interaction_genes[idx[0]])

		set_length = len(list(set(final_genes)))

	return final_genes


# get a random sample of genes from interaction_genes that DO NOT have any
# interactions with one another
# num_genes: number of random genes to select
# total_gene_list: list of genes within the dataset being used
# interaction_genes: list of valid genes from an interaction list
# interaction_list: should be two columns: interacting A and interacting B...
#					a list of interacting genes (pairwise)
# NOTE: it is assumed that interaction_genes are only genes contained inside of
# total_gene_list, perform a check before passing it to this function
def create_random_subset_from_NON_interactions(num_genes, total_gene_list, \
											interaction_genes, interaction_list):
	final_genes = []

	# start a list with a single random gene
	gene_indexes = np.random.randint(0, len(interaction_genes), 1)
	final_genes = [interaction_genes[gene_indexes[0]]]

	# loop through and add a random gene if it is not already in the list
	# and if it is does not have any interactions
	while len(final_genes) != num_genes:
		interact = 0
		gene_indexes = np.random.randint(0, len(interaction_genes), 1)
		rand_gene = interaction_genes[gene_indexes[0]]

		if rand_gene in final_genes:
			continue

		locs = np.where(interaction_list == rand_gene)

		for i in range(locs[0].shape[0]):
			cand_gene = interaction_list[locs[0][i], 1 - locs[1][i]]
			if cand_gene in final_genes:
				interact = 1
				break

		if interact:
			continue

		final_genes.append(rand_gene)

	return final_genes


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

		for i in range(1, set_size):
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

	# only use human genes
	biogrid = biogrid[biogrid['ORGANISM_A_ID'] == 9606]
	biogrid = biogrid[biogrid['ORGANISM_B_ID'] == 9606]

	# only interested in offical symbols for A and B and when genes are not equal
	interactions = biogrid[biogrid['OFFICIAL_SYMBOL_A'] != biogrid['OFFICIAL_SYMBOL_B']]
	interactions = interactions[['OFFICIAL_SYMBOL_A', 'OFFICIAL_SYMBOL_B']]
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
					print(gene)
					if g not in interacted_subsets[s]:
						interacted_subsets[s].append(g)
						print(g)

	return interacted_subsets


# USAGE:
# - return a list of subsets only containing genes with NO interactions within the subset
# PARAMS:
# subsets - dictionary of subsets similar to whats returned by read_subset_file
# interaction_file - file containing gene interactions

# NOTE: the interaction file should not contain any metadata in the top (i.e. a readme)...
# it should just contain the column headers followed by the interaction data.
# column format is expected to follow the biogrid interaction file 
def convert_subset_to_non_interactions(subsets, interaction_file):
	# read biogrid file into a dataframe, then extract official symbols
	biogrid = pd.read_csv(interaction_file, sep='\t')

	# only use human genes
	biogrid = biogrid[biogrid['ORGANISM_A_ID'] == 9606]
	biogrid = biogrid[biogrid['ORGANISM_B_ID'] == 9606]

	# only interested in offical symbols for A and B and when genes are not equal
	interactions = biogrid[biogrid['OFFICIAL_SYMBOL_A'] != biogrid['OFFICIAL_SYMBOL_B']]
	interactions = interactions[['OFFICIAL_SYMBOL_A', 'OFFICIAL_SYMBOL_B']]
	interactions = interactions.values
	interactions = interactions.astype(np.str)

	non_interacted_subsets = {}
	for s in subsets:
		non_interacted_subsets[s] = []
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

			int_gs = []
			for gene in total_interacting_genes:
				if gene in subsets[s]:
					int_gs.append(gene)

			if len(int_gs) == 0:
				non_interacted_subsets[s].append(g)


	return non_interacted_subsets


# write a subset list dictionary to a file
def write_subsets_to_file(subsets, file):
	with open(file, 'w') as f:
		for s in subsets:
			f.write(s + '\t')
			for g in subsets[s]:
				f.write(g + '\t')
			f.write('\n')
