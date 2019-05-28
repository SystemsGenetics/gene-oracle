import argparse
import dataframe_helper
import json
import numpy as np
import pandas as pd
import random
import sklearn.dummy
import sklearn.model_selection
import sklearn.preprocessing
import sys



def load_gene_sets(filename):
	# load file into list
	lines = [line.strip() for line in open(filename, "r")]
	lines = [line.split("\t") for line in lines]

	# map each gene set into a tuple of the name and genes in the set
	gene_sets = [(line[0], line[1:]) for line in lines]

	return gene_sets



def evaluate(data, labels, clf, genes, cv=5):
	# extract dataset
	X = data[genes]

	# normalize dataset
	X = sklearn.preprocessing.MaxAbsScaler().fit_transform(X)

	# evaluate gene set
	scores = sklearn.model_selection.cross_val_score(clf, X, y=labels, cv=cv)

	return list(scores)



def evaluate_set(data, labels, clf, name, genes, cv=5, outfile=None):
	# evaluate gene set
	scores = evaluate(data, labels, clf, genes, cv=cv)

	# print results
	line = "\t".join([name] + ["%.3f" % (score) for score in scores])

	print(line)

	# write results to output file
	if outfile:
		outfile.write(line + "\n")



def evaluate_random(data, labels, clf, n_genes, cv=5, n_iters=100, outfile=None):
	# evaluate n_iters random sets
	scores = []

	for i in range(n_iters):
		# generate random gene set
		genes = random.sample(list(data.columns), n_genes)

		# evaluate gene set
		scores += evaluate(data, labels, clf, genes, cv=cv)

	# print results
	name = "random-%d" % n_genes
	line = "\t".join([name] + ["%.3f" % (score) for score in scores])

	print(line)

	# write results to output file
	if outfile:
		outfile.write(line + "\n")



if __name__ == "__main__":
	# parse command-line arguments
	parser = argparse.ArgumentParser(description="Evaluate classification potential of gene sets")
	parser.add_argument("--dataset", help="input dataset (samples x genes)", required=True)
	parser.add_argument("--labels", help="list of sample labels", required=True)
	parser.add_argument("--model_config", help="json file containing network specifications", required=True)
	parser.add_argument("--outfile", help="output file to save results")
	parser.add_argument("--gene_sets", help="list of gene sets to evaluate (GMT/GCT format)")
	parser.add_argument("--full", help="Evaluate the set of all genes in the dataset", action="store_true")
	parser.add_argument("--random", help="Evaluate random gene sets", action="store_true")
	parser.add_argument("--random_range", help="range of random gene sizes to evaluate", nargs=2, type=int)
	parser.add_argument("--random_iters", help="number of iterations to perform for random classification", type=int, default=100)
	parser.add_argument("--num_folds", help="number of folds for k-fold cross validation", type=int, default=5)

	args = parser.parse_args()

	# load input data
	print("loading input dataset...")

	df = dataframe_helper.load(args.dataset)
	df_samples = df.index
	df_genes = df.columns
	labels = pd.read_csv(args.labels, sep="\t", header=None, index_col=0)

	print("loaded input dataset (%s genes, %s samples)" % (df.shape[1], df.shape[0]))

	# initialize classifier
	print("initializing classifier...")

	model_config = json.load(open(args.model_config))

	# TODO: initialize MLP
	clf = sklearn.dummy.DummyClassifier()

	# load subset file if it was provided
	if args.gene_sets != None:
		print("loading gene sets...")

		gene_sets = load_gene_sets(args.gene_sets)

		print("loaded %d gene sets" % (len(gene_sets)))

		# remove genes which do not exist in the dataset
		genes = list(set(sum([genes for (name, genes) in gene_sets], [])))
		missing_genes = [g for g in genes if g not in df_genes]

		gene_sets = [(name, [g for g in genes if g in df_genes]) for (name, genes) in gene_sets]

		print("%d / %d (%0.1f%%) genes from gene sets were not found in the input dataset" % (
			len(missing_genes),
			len(genes),
			len(missing_genes) / len(genes) * 100))
	else:
		gene_sets = []

	# include the set of all genes if specified
	if args.full:
		gene_sets.append(("FULL", df_genes))

	# initialize list of random set sizes
	if args.random:
		# determine random set sizes from range
		if args.random_range != None:
			print("initializing random set sizes from range...")
			random_sets = range(args.random_range[0], args.random_range[1] + 1)

		# determine random set sizes from gene sets
		elif args.gene_sets != None:
			print("initializing random set sizes from gene sets...")
			random_sets = sorted(set([len(genes) for (name, genes) in gene_sets]))

		# print error and exit
		else:
			print("error: --gene_sets or --random_range must be provided to determine random set sizes")
			sys.exit(1)
	else:
		random_sets = []

	print("evaluating gene sets...")

	# write header to output file
	if args.outfile:
		outfile = open(args.outfile, "w")
		outfile.write("\t".join(["Name"] + ["%d" % (i) for i in range(args.num_folds)]) + "\n")

	# evaluate input gene sets
	for (name, genes) in gene_sets:
		evaluate_set(df, labels[1], clf, name, genes, cv=args.num_folds, outfile=outfile)

	# evaluate random gene sets
	for n_genes in random_sets:
		evaluate_random(df, labels[1], clf, n_genes, cv=args.num_folds, n_iters=args.random_iters, outfile=outfile)
