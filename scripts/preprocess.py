"""
This script performs various pre-processing tasks on gene expression matrix (GEM)
files, including (1) conversion between plaintext and numpy format and (2) transposition.
"""
import argparse
import numpy as np
import pandas as pd
import sys

import utils



if __name__ == "__main__":
	# parse command-line arguments
	parser = argparse.ArgumentParser(description="Preprocess a GEM to be compatible with gene-oracle")
	parser.add_argument(dest="infile", help="Input file")
	parser.add_argument(dest="outfile", help="Output file")
	parser.add_argument("--transpose", help="Transpose the input file", action="store_true")

	args = parser.parse_args()

	# load dataframe from input format
	print("loading %s..." % args.infile)

	df = utils.load_dataframe(args.infile)

	# transpose dataframe if specified
	if args.transpose:
		print("transposing matrix...")
		df = df.T

	# save dataframe in output format
	print("saving %s..." % args.outfile)

	utils.save_dataframe(args.outfile, df)
