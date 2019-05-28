import argparse
import json
import numpy as np
import sys



def load_gene_sets(filename):
	# load file into list
	lines = [line.strip() for line in open(filename, "r")]
	lines = [line.split("\t") for line in lines]

	# map each gene set into a tuple of the name and genes in the set
	gene_sets = [(line[0], line[1:]) for line in lines]

	return gene_sets



def evaluate_full():
	pass



def evaluate_random():
	pass



def evaluate_set():
	pass



if __name__ == "__main__":
	# parse command-line arguments
	parser = argparse.ArgumentParser(description="Run classification on specified dataset, subset of genes, or a random set")
	parser.add_argument("--data", help="numpy array of input dataset", required=True)
	parser.add_argument("--rownames", help="list of genes (rows in dataset)", required=True)
	parser.add_argument("--colnames", help="list of samples (columns in dataset)", required=True)
	parser.add_argument("--model_config", help="json file containing network specifications", required=True)
	parser.add_argument("--outfile", help="output file to save results", required=True)
	parser.add_argument("--gene_sets", help="list of gene sets to evaluate (GMT/GCT format)")
	parser.add_argument("--random", help="Evaluate random gene sets", action="store_true")
	parser.add_argument("--random_range", help="range of random gene sizes to evaluate", nargs=2, type=int)
	parser.add_argument("--random_iters", help="number of iterations to perform for random classification", type=int, default=100)
	parser.add_argument("--num_folds", help="number of folds for k-fold cross validation", type=int, default=5)

	args = parser.parse_args()

	# TODO: validate args
	
	# load input data
	print("loading input dataset...")

	data = np.load(args.data)
	rownames = np.loadtxt(args.rownames, dtype=str)
	colnames = np.loadtxt(args.colnames, dtype=str, delimiter="\t")

	# validate input data
	if data.shape[0] != len(rownames):
		print("error: data and rownames do no match")
		sys.exit(1)

	if data.shape[1] != len(colnames):
		print("error: data and colnames do no match")
		sys.exit(1)

	if colnames.shape[1] != 2:
		print("error: colnames should have two columns, one for column names and one for labels")
		sys.exit(1)

	print("loaded input dataset (%s genes, %s samples)" % (len(rownames), len(colnames)))

	# load model config
	print("loading model config...")

	model_config = json.load(open(args.model_config))

	# load subset file if it was provided
	if args.gene_sets != None:
		print("loading gene sets...")

		gene_sets = load_gene_sets(args.gene_sets)

		print("loaded %d gene sets" % (len(gene_sets)))

		# remove genes which do not exist in the dataset
		genes = list(set(sum([genes for (name, genes) in gene_sets], [])))
		missing_genes = [g for g in genes if g not in rownames]

		gene_sets = [(name, [g for g in genes if g in rownames]) for (name, genes) in gene_sets]

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
		print("evaluating global gene set...")

		evaluate_full()

	# evaluate input gene sets
	print("evaluating input gene sets...")

	for (name, genes) in gene_sets:
		print("  evaluating %s (%d genes)..." % (name, len(genes)))

		evaluate_set()

	# evaluate random gene sets
	print("evaluating random gene sets...")

	for size in random_sets:
		print("  evaluating random %d..." % (size))

		evaluate_random()
