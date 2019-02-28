import argparse
import pandas as pd
import scipy.stats

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
	lines = [line.split(",") for line in lines]
	subset_dict = {line[0]: line[1:] for line in lines}

	# evaluate each curated gene set
	print("%s\t%s\t%s" % ("Name", "t", "p"))

	for subset_name in df_subset.index:
		subset_size = len(subset_dict[subset_name])

		t, p = scipy.stats.ttest_ind(df_subset.loc[subset_name], df_random.loc[subset_size], equal_var=False)

		print("%s\t%f\t%f" % (subset_name, t, p))
