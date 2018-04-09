'''

delta_accs_box.py

	This script allows for plotting of the differences between the classification accuracies of random genes ve the
	classification accuracy of hallmark genes

	The program expects the following parameter
	--rand_dir : dir containing random accuracies as a string
	--sub_dir : dir containing subset accuracies as a string


	example: python delta_accs_box.py --rand_accs "../logs/panTCGA_hallmark_random.log" --sub_accs "../logs/panTCGA_hallmark_results.log" --sub_count "../subsets/hallmark_gene_counts.json" --out "test.txt"
	

'''
import sys, argparse
import matplotlib.pyplot as plt
import json
import numpy as np
import seaborn as sns
import operator as op
import pandas as pd

# read in the accuracy files
def read_file_avg(file):
	avg_accs = {}
	with open(file, 'r') as f:
		next(f) # skip the first line which has string header info
		for line in f:
			line = line.split('\t')
			avg_accs[line[0]] = float(line[1])

	return avg_accs

def read_file(file):
	accs = {}
	with open(file, 'r') as f:
		next(f) # skip the first line which has string header info
		for line in f:
			line = line.split('\t')
			info = []
			info.append(float(line[1]) + float(line[2]))# the average plus the standarddev
			info.append(float(line[1]) - float(line[2]))# the average minus the standarddev
			info.append(float(line[1]))#avg
			info.append(float(line[3]))#max
			info.append(float(line[4]))#min
			accs[line[0]] = info

	return accs


# plot the delta accuracies
def plot(d_accs, gene_counts, out):

	fig, ax = plt.subplots()
	fig.set_size_inches(50,50)
	sorted_keys, sorted_vals = zip(*sorted(d_accs.items(),key=op.itemgetter(1)))
	ax = sns.boxplot(data = sorted_vals,orient='h')
	ax.set_yticklabels(delta_accs.keys())
	plt.show()
	

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Visualize delta accuracy between random genes and subsets')
	parser.add_argument('--rand_accs', help='file random accuracies', type=str, required=True)
	parser.add_argument('--sub_accs', help='file of subset accuracies', type=str, required=True)
	parser.add_argument('--sub_count', help='json with count of genes in subset', type=str, required=True)
	parser.add_argument('--out', help='file to save to', type=str, required=False)
	# additional args?
	args = parser.parse_args()
	
	# read two accuracy files
	rand_accs = read_file(args.rand_accs)
	sub_accs = read_file_avg(args.sub_accs)

	# read a json file containing the count of genes for each subset
	with open(args.sub_count, 'r') as f:
		gene_count_dict = json.load(f)


	print(rand_accs)
	print(sub_accs)

	# get the delta accs
	gene_counts = []
	delta_accs = {}
	for s in sorted(gene_count_dict.keys()):
		delta_accs[s] = []

		for val in rand_accs[str(gene_count_dict[s])]:

			delta = round(sub_accs[s] - val, 3)

			if delta < 0:
				delta = 0.0

			delta_accs[s].append(delta)

			gene_counts.append(gene_count_dict[s])

	for key, value in delta_accs.items():
		print(str(key) + " , " + str(value) )


	plot(delta_accs, gene_counts, args.out)