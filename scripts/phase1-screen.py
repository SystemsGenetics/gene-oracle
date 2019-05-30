"""
This script uses Welch's t-test to determine whether a gene set exhibits
significantly higher classification accuracy over equivalent random sets. The
script takes three inputs -- a list of gene sets and their accuracies, a list of
random sets and their accuracies, and a list of gene sets and their member
genes -- and produces a list of t-test results for each gene set. The lower the
p-value, the more likely it is that the corresponding gene set performs
significantly better over random.
"""
import argparse
import pandas as pd
import scipy.stats

import utils



def evaluate_from_samples(df_subset, df_random, name, n_genes):
	return scipy.stats.ttest_ind(df_subset.loc[name], df_random.loc[n_genes], equal_var=False)



def evaluate_from_stats(df_subset, n_subset, df_random, n_random, name, n_genes):
	return scipy.stats.ttest_ind_from_stats( \
		df_subset.loc[name, "Average"], df_subset.loc[name, "Std_Dev"], n_subset, \
		df_random.loc[n_genes, "Average"], df_random.loc[n_genes, "Std_Dev"], n_random, \
		equal_var=False)



if __name__ == "__main__":
	# parse command-line arguments
	parser = argparse.ArgumentParser(description="Select gene sets which perform significantly better than equivalent random sets.")
	parser.add_argument("--random", help="list of accuracies for random gene sets", required=True)
	parser.add_argument("--subset", help="list of accuracies for curated gene sets", required=True)
	parser.add_argument("--gene_sets", help="list of gene sets (GMT/GCT)")

	args = parser.parse_args()

	# load input files
	df_random = pd.read_csv(args.random, sep="\t", header=None, index_col=0)
	df_subset = pd.read_csv(args.subset, sep="\t", header=None, index_col=0)
	gene_sets = utils.load_gene_sets(args.gene_sets)

	# evaluate each curated gene set
	print("%-80s %s" % ("Name", "p"))

	for name, genes in gene_sets:
		t, p = evaluate_from_samples(df_subset, df_random, name, len(genes))

		print("%-80s %0.12f" % (name, p))
