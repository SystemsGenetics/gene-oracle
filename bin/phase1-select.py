#!/usr/bin/env python

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
import operator
import pandas as pd
import scipy.stats

import utils



if __name__ == "__main__":
	# parse command-line arguments
	parser = argparse.ArgumentParser(description="Select gene sets which perform significantly better than equivalent random sets.")
	parser.add_argument("--scores", help="list of scores for curated and random gene sets", required=True)
	parser.add_argument("--gene-sets", help="list of curated gene sets")
	parser.add_argument("--threshold", help="maximum p-value required for a gene set to be selected", type=float, default=1)
	parser.add_argument("--n-sets", help="maximum number of gene sets that can be selected", type=int, default=-1)
	parser.add_argument("--outfile", help="output file to save results", required=True)

	args = parser.parse_args()

	# load input files
	scores = pd.read_csv(args.scores, sep="\t", index_col=0)
	gene_sets = utils.load_gene_sets(args.gene_sets)

	# evaluate each curated gene set
	results = []

	for name, genes in gene_sets:
		# perform t-test between gene set score and background score
		index_fg = name
		index_bg = str(len(genes))

		t, p = scipy.stats.ttest_ind_from_stats( \
			scores.loc[index_fg, "mu"], scores.loc[index_fg, "sigma"], scores.loc[index_fg, "n"], \
			scores.loc[index_bg, "mu"], scores.loc[index_bg, "sigma"], scores.loc[index_bg, "n"], \
			equal_var=False)

		results.append((name, genes, p))

		# print result
		print("%-40s %0.12f" % (name, p))

	# select gene sets
	results.sort(key=operator.itemgetter(2))
	results = results[0:args.n_sets]
	results = [(name, genes) for (name, genes, p) in results if p <= args.threshold]

	# save selected gene sets to output file
	outfile = open(args.outfile, "w")

	for (name, genes) in results:
		outfile.write("\t".join([name] + genes) + "\n")
