import argparse
import dataframe_helper
import json
import numpy as np
import pandas as pd
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



def evaluate_random():
	pass



def evaluate_set(data, labels, clf, name, genes, cv=5, outfile=None):
	# extract dataset
	X = data[genes]

	# normalize dataset
	X = sklearn.preprocessing.MaxAbsScaler().fit_transform(X)

	# evaluate gene set
	scores = sklearn.model_selection.cross_val_score(clf, X, y=labels, cv=cv)

	# print results
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
	parser.add_argument("--random", help="Evaluate random gene sets", action="store_true")
	parser.add_argument("--random_range", help="range of random gene sizes to evaluate", nargs=2, type=int)
	parser.add_argument("--random_iters", help="number of iterations to perform for random classification", type=int, default=100)
	parser.add_argument("--num_folds", help="number of folds for k-fold cross validation", type=int, default=5)

	args = parser.parse_args()

	# TODO: validate args
	
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

	# initialize list of random set sizes
	if args.random:
		if args.random_range != None:
			print("initializing random set sizes from range...")
			random_sets = range(args.random_range[0], args.random_range[1] + 1)
		elif args.gene_sets != None:
			print("initializing random set sizes from gene sets...")
			random_sets = sorted(set([len(genes) for (name, genes) in gene_sets]))
		else:
			print("error: --random was specified without any way to determine random set sizes")
			sys.exit(1)
	else:
		random_sets = []

	# evaluate the set of all genes if no gene sets or random sets were provided
	if args.gene_sets == None and not args.random:
		gene_sets.append(("FULL", df_genes))

	# write header to output file
	if args.outfile:
		outfile = open(args.outfile, "w")
		outfile.write("\t".join(["Name"] + ["%d" % (i) for i in range(args.num_folds)]) + "\n")

	# evaluate input gene sets
	print("evaluating input gene sets...")

	for (name, genes) in gene_sets:
		evaluate_set(df, labels[1], clf, name, genes, cv=args.num_folds, outfile=outfile)

	# evaluate random gene sets
	print("evaluating random gene sets...")

	for size in random_sets:
		print("  evaluating random %d..." % (size))

		evaluate_random()
