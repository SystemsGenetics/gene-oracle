#
# This file takes in a text file containing genetic subsets in
# symbolic form and outputs the same subsets in ENSG format
#

import numpy as np
import argparse
from utils import read_subset_file


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Convert subsets to ENSG format')
	parser.add_argument('--subset_list', help='gmt/gct/txt file containing subsets', type=str, required=True)
	parser.add_argument('--gene_map', help='file containing gene mapping', type=str, required=True)
	parser.add_argument('--out_file', help='file to write new list to', type=str, required=True)
	args = parser.parse_args()

	# read the original subset
	orig_subs = read_subset_file(args.subset_list)

	# read in gene mapping file
	gm = np.loadtxt(args.gene_map, dtype=np.str)
	
	new_gm = []
	# format the gene mapping
	for i in range(gm.shape[0]):
		new_gm.append(gm[i].split(','))
	gm = np.asarray(new_gm)

	# instantiate the new dictionary for ensg format
	ensg_subs = {}

	tot_genes = []
	missing_genes = []

	# iterate through each subset
	for s in orig_subs:
		ensg_subs[s] = []
		missing = 0

		# look for each gene in the gene map
		for g in orig_subs[s]:
			if g not in tot_genes:
				tot_genes.append(g)

			idxs = np.where(gm[:,2] == g)
			if idxs[0].shape[0]:
				ensg_subs[s].append(gm[idxs[0][0], 0])
			else:
				missing += 1
				if g not in missing_genes:
					missing_genes.append(g)
		print(str(s) + ' missing ' + str(missing) + ' out of ' + \
										str(len(orig_subs[s]))+ ' genes')

	# write the new output file
	with open(args.out_file, 'w') as f:
		for s in ensg_subs:
			f.write(str(s) + ',')
			for g in ensg_subs[s][:-1]:
				f.write(str(g) + ',')
			f.write(str(ensg_subs[s][-1]) + '\n')

	# done!
	print('done mapping')
	print('missing ' + str(len(missing_genes)) + '/' + str(len(tot_genes)) + ' or ' + \
		str(int((float(len(missing_genes)) / len(tot_genes)) * 100.0)) + '% of genes')





