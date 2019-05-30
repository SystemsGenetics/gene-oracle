import argparse
import utils
import numpy as np
import pandas as pd
import random
import sklearn.datasets
import sys



if __name__ == "__main__":
	# parse command-line arguments
	parser = argparse.ArgumentParser(description="Create a synthetic classification dataset")
	parser.add_argument("--n-samples", help="number of samples", type=int, default=100)
	parser.add_argument("--n-genes", help="number of genes", type=int, default=20)
	parser.add_argument("--n-classes", help="number of classes", type=int, default=2)
	parser.add_argument("--n-sets", help="number of gene sets", type=int, default=10)

	args = parser.parse_args()

	# create synthetic dataset
	n_informative = args.n_genes // 10
	n_redundant = args.n_genes - n_informative

	X, y = sklearn.datasets.make_classification(args.n_samples, args.n_genes, n_informative=n_informative, n_redundant=n_redundant, n_classes=args.n_classes)

	# initialize dataframe
	X = pd.DataFrame(X)
	X_genes = [str(i) for i in X.columns]
	y = pd.DataFrame(y)

	# create synthetic gene sets
	gene_sets = []

	for i in range(args.n_sets):
		n_genes = random.randint(3, args.n_genes)
		genes = random.sample(X_genes, n_genes)

		gene_sets.append(["gene-set-%d" % i] + genes)

	# save dataset to file
	utils.save_dataframe("example_data.txt", X)

	# save labels to file
	y.to_csv("example_labels.txt", sep="\t", header=None)

	# save gene sets to file
	f = open("example_genesets.txt", "w")
	f.write("\n".join(["\t".join(gene_set) for gene_set in gene_sets]))
