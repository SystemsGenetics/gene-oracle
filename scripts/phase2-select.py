"""
This script takes the subset scores for a gene set and measures the saliency of
each individual gene by how frequently it occurs in all subsets that were
evaluated. Genes with a higher "aggregate frequency" are selected as "candidate"
genes, while the other genes are labeled as "non-candidate" genes.
"""
import argparse
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import sklearn.mixture

import utils



def load_subsets(logdir, name, n_genes):
	# load all subset scores for a gene set
	subsets = []

	for k in range(1, n_genes + 1):
		logfile = open("%s/%s_scores_%03d.txt" % (logdir, name, k), "r")
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



def compute_scores(freq_matrix):
	return freq_matrix.sum(axis=0)



def compute_threshold(genes, freq_matrix):
	# compute aggregate frequency of each gene
	scores = freq_matrix.sum(axis=0)

	# fit a Gaussian mixture model to the gene scores
	X = scores.reshape(-1, 1)

	gmm = sklearn.mixture.GaussianMixture(n_components=2)
	gmm.fit(X)

	# compute the intersection between the two modes
	m1 = gmm.means_[0, 0]
	m2 = gmm.means_[1, 0]
	s1 = gmm.covariances_[0, 0, 0]
	s2 = gmm.covariances_[1, 0, 0]

	num = m2*s1**2 - m1*s2**2
	delta = s1 * s2 * math.sqrt((m1 - m2)**2 + 2 * (s1**2 - s2**2) * math.log(s1/s2))
	denom = s1**2 - s2**2

	m = (m1 + m2) / 2
	c1 = (num + delta) / denom
	c2 = (num - delta) / denom

	if abs(c1 - m) < abs(c2 - m):
		threshold = c1
	else:
		threshold = c2

	return threshold, scores



if __name__ == "__main__":
	# parse command-line arguments
	parser = argparse.ArgumentParser(description="Identify candidate / non-candidate genes in a gene set")
	parser.add_argument("--gene-sets", help="list of curated gene sets", required=True)
	parser.add_argument("--logdir", help="directory where logs are stored", required=True)
	parser.add_argument("--threshold", help="manual threshold based on percentile (0-100)", type=float)
	parser.add_argument("--visualize", help="visualize frequency heatmap and candidate threshold", action="store_true")
	parser.add_argument("--outfile", help="output file to save results", required=True)

	args = parser.parse_args()

	# load gene sets
	gene_sets = utils.load_gene_sets(args.gene_sets)

	# initialize output file
	outfile = open(args.outfile, "w")

	# select candidate genes for each gene set
	for name, genes in gene_sets:
		# compute frequency matrix
		subsets = load_subsets(args.logdir, name, len(genes))
		freq_matrix = compute_frequency_matrix(genes, subsets)

		# plot heatmap of frequency matrix
		if args.visualize:
			sns.heatmap(freq_matrix, xticklabels=genes, yticklabels=range(1, len(genes) + 1))
			plt.title(name)
			plt.savefig("%s-frequency-heatmap.png" % (name))
			plt.close()

		# compute aggregate frequency of each gene
		scores = compute_scores(freq_matrix)

		# use a percentile threshold if specified
		if args.threshold != None:
			threshold = np.percentile(scores, args.threshold)

		# otherwise compute threshold automatically
		else:
			threshold = compute_threshold(genes, scores)

		# select candidate genes
		candidate_genes = [gene for i, gene in enumerate(genes) if scores[i] > threshold]

		# plot distribution of gene scores
		if args.visualize:
			sns.distplot(scores)
			ymin, ymax = plt.gca().get_ylim()
			y = [ymin, ymax / 2]
			plt.plot([threshold, threshold], y, "r")
			plt.title(name)
			plt.savefig("%s-candidate-threshold.png" % (name))
			plt.close()

		# save results to output file
		outfile.write("\t".join([name] + candidate_genes) + "\n")
