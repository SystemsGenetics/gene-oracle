"""
This script uses Welch's t-test to determine whether a gene set exhibits
significantly higher classification accuracy over equivalent random sets. The
script takes three inputs -- a list of curated sets and their scores, a list of
random sets and their scores, and a list of gene sets and their member
genes -- and produces a list of t-test results for each gene set. The lower the
p-value, the more likely it is that the corresponding gene set performs
significantly better over random.
"""
import argparse
import pandas as pd
import scipy.stats

import utils



def evaluate_from_samples(data_fg, index_fg, data_bg, index_bg):
	return scipy.stats.ttest_ind( \
		data_fg.loc[index_fg], \
		data_bg.loc[index_bg], \
		equal_var=False)



def evaluate_from_stats(data_fg, index_fg, data_bg, index_bg):
	return scipy.stats.ttest_ind_from_stats( \
		data_fg.loc[index_fg, "mu"], data_fg.loc[index_fg, "sigma"], data_fg.loc[index_fg, "n"], \
		data_bg.loc[index_bg, "mu"], data_bg.loc[index_bg, "sigma"], data_bg.loc[index_bg, "n"], \
		equal_var=False)



if __name__ == "__main__":
	# parse command-line arguments
	parser = argparse.ArgumentParser(description="Select gene sets which perform significantly better than equivalent random sets.")
	parser.add_argument("--scores-fg", help="list of scores for curated gene sets (foreground)", required=True)
	parser.add_argument("--scores-bg", help="list of scores for random gene sets (background)", required=True)
	parser.add_argument("--gene-sets", help="list of gene sets")

	args = parser.parse_args()

	# load input files
	data_fg = pd.read_csv(args.scores_fg, sep="\t", header=None, index_col=0)
	data_bg = pd.read_csv(args.scores_bg, sep="\t", header=None, index_col=0)
	gene_sets = utils.load_gene_sets(args.gene_sets)

	# evaluate each curated gene set
	print("%-80s %s" % ("Name", "p"))

	for name, genes in gene_sets:
		t, p = evaluate_from_samples(data_fg, name, data_bg, len(genes))

		print("%-80s %0.12f" % (name, p))
