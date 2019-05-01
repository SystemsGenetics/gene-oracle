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



def evaluate_samples(df_subset, df_random, subset_name, subset_size):
	return scipy.stats.ttest_ind(df_subset.loc[subset_name], df_random.loc[subset_size], equal_var=False)



def evaluate_stats(df_subset, df_random, subset_name, subset_size):
	return scipy.stats.ttest_ind_from_stats( \
		df_subset.loc[subset_name, "Average"], df_subset.loc[subset_name, "Std_Dev"], 10, \
		df_random.loc[subset_size, "Average"], df_random.loc[subset_size, "Std_Dev"], 500, \
		equal_var=False)



if __name__ == "__main__":
	# parse command-line arguments
	parser = argparse.ArgumentParser(description="Select gene sets which perform significantly better than equivalent random sets.")
	parser.add_argument("--random", help="Log file of accuracies for random sets", required=True)
	parser.add_argument("--subset", help="Log file of accuracies for specific gene sets", required=True)
	parser.add_argument("--dict", help="Text file of gene lists", required=True)

	args = parser.parse_args()

	# load input files
	df_random = pd.read_csv(args.random, sep="\t", index_col=0)
	df_subset = pd.read_csv(args.subset, sep="\t", index_col=0)
	
	# load gene set dictionary
	lines = [line.rstrip() for line in open(args.dict, "r")]
	lines = [line.split("\t") for line in lines]
	subset_dict = {line[0]: line[1:] for line in lines}

	# evaluate each curated gene set
	print("%-80s %s" % ("Name", "p"))

	for subset_name in df_subset.index:
		subset_size = len(subset_dict[subset_name])

		t, p = evaluate_samples(df_subset, df_random, subset_name, subset_size)

		print("%-80s %0.12f" % (subset_name, p))
