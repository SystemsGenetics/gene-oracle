"""
This script takes the subset scores for a gene set and measures the saliency of
each individual gene by how frequently it occurs in all subsets that were
evaluated. Genes with a higher "aggregate frequency" are selected as "candidate"
genes, while the other genes are labeled as "non-candidate" genes.
"""
import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

import utils



def load_subsets(logdir, n_genes):
	# load all subset scores for a gene set
	subsets = []

	for k in range(1, n_genes + 1):
		logfile = open("%s/scores_%03d.txt" % (logdir, k), "r")
		lines = [line.strip() for line in logfile]
		lines = [line.split("\t") for line in lines]

		subsets += [(line[0].split(","), float(line[1])) for line in lines]

	return subsets



def compute_frequency_matrix(genes, subsets):
	# initialize frequency matrix
	n_genes = len(genes)
	freq_matrix = np.zeros((n_genes, n_genes))

	# initialize dictionary of gene indices for quick access
	gene_dict = {gene: genes.index(gene) for gene in genes}

	# populate frequency matrix from subsets
	for subset_genes, subset_score in subsets:
		k = len(subset_genes)

		for gene in subset_genes:
			freq_matrix[k - 1, gene_dict[gene]] += 1

	# normalize freqency matrix by the number of genes in each iteration
	freq_matrix /= freq_matrix.sum(axis=1, keepdims=True)

	return freq_matrix



def select_candidate_genes(genes, freq_matrix):
	# compute aggregate frequency of each gene
	frequency_sums = freq_matrix.sum(axis=0)

	# plot distribution of aggregate frequencies
	sns.distplot(frequency_sums)
	plt.show()



if __name__ == "__main__":
	# parse command-line arguments
	parser = argparse.ArgumentParser(description="Identify candidate / non-candidate genes in a gene set")
	parser.add_argument("--gene-sets", help="list of curated gene sets", required=True)
	parser.add_argument("--logdir", help="directory where logs are stored", required=True)

	args = parser.parse_args()

	# load gene sets
	gene_sets = utils.load_gene_sets(args.gene_sets)

	# select candidate genes for each gene set
	for name, genes in gene_sets:
		print(name)

		# compute frequency matrix
		logdir = "%s/%s" % (args.logdir, name)
		subsets = load_subsets(logdir, len(genes))
		freq_matrix = compute_frequency_matrix(genes, subsets)

		# plot heatmap of frequency matrix
		sns.heatmap(freq_matrix, xticklabels=genes, yticklabels=range(1, len(genes) + 1))
		plt.title(name)
		plt.show()

		select_candidate_genes(genes, freq_matrix)
