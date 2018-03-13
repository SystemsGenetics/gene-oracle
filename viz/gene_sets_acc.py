'''

gene_sets_acc.py

	This script allows for plotting of a log directory that contains accuracies subsets of
	a particular set of genes. Additionally the program accepts a directory containing
	random gene sets so you can plot multiple runs at once

	The program expects a file that contains the following format:

	[list of genes] accuracy
	...

'''



import sys, argparse
import matplotlib.pyplot as plt
import json
import numpy as np
import seaborn as sns
import pandas as pd
from scipy import stats
import os

# read in the accuracy files to dictionary
def read_sub_dir(s_dir):
	accs = {}
	for sub_f in os.listdir(s_dir):
		sub_f_split = sub_f.split('_')
		idx = int([s for s in sub_f_split if s.isdigit()][0])
		accs[idx] = []
		with open(os.path.join(s_dir, sub_f), 'r') as f:
			for line in f:
				line = line.split('\t')
				accs[idx].append(float(line[1]))
	return accs


# read in sub dir accuracies from a top level dir
def read_top_dir(t_dir):
	accs_list = []

	for sub_d in os.listdir(t_dir):
		accs_list.append(read_sub_dir(os.path.join(t_dir,sub_d)))

	return accs_list


# plot the delta accuracies
def plot(rand_accs, sub_accs, rand_std_es, sub_stes, avg_rands):
	
	avgs = []
	for s in sorted(sub_accs.keys()):
		avgs.append(sum(sub_accs[s]) / len(sub_accs[s]))

	num = sorted(sub_accs.keys())

	(_, caps, _) = plt.errorbar(num, avgs, sub_stes, fmt='ok', lw=1, markersize=5, markerfacecolor='b', capsize=8)

	for cap in caps:
	    cap.set_color('b')
	    cap.set_markeredgewidth(.5)

	if avg_rands:
		avgs = []
		stes = []
		for k in rand_accs[0]:
			concat = []
			for r in rand_accs:
				concat += r[k]
			avgs.append(sum(concat) / len(concat))
			stes.append(stats.sem(concat))

		(_, caps, _) = plt.errorbar(num, avgs, stes, fmt='ok', lw=1, markersize=5, markerfacecolor='r', capsize=8)

		for cap in caps:
		    cap.set_color('r')
		    cap.set_markeredgewidth(.5)

	else:
		colors = ["r","b","g","c","m","y","k","b","m","r"]

		for r, s, c in zip(rand_accs, rand_std_es, colors):
			avgs = []
			for a in sorted(r.keys()):
				avgs.append(sum(r[a]) / len(r[a]))

			(_, caps, _) = plt.errorbar(num, avgs, s, fmt='ok', lw=1, markersize=5, markerfacecolor=c, capsize=8)

			for cap in caps:
			    cap.set_color(c)
			    cap.set_markeredgewidth(.5)

	plt.title('Gene Set Accuracies')
	plt.xlabel('Num Genes in Set')
	plt.ylabel('Accuracy')
	plt.xticks(num)

	ax = plt.axes()
	ax.xaxis.grid(linestyle='--')

	plt.show()


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Visualize delta accuracy between random genes and subsets')
	parser.add_argument('--rand_dir', help='dir containing random accuracies', type=str, required=True)
	parser.add_argument('--sub_dir', help='dir containing subset accuracies', type=str, required=True)
	parser.add_argument('--avg_rands', help='average the random runs', action='store_true', required=False)
	# additional args?
	args = parser.parse_args()

	# read two accuracy files
	rand_accs = read_top_dir(args.rand_dir)
	sub_accs = read_sub_dir(args.sub_dir)

	# calculate standard error
	rand_std_es = []
	for r in rand_accs:
		std_es = []
		for s in r:
			std_es.append(stats.sem(r[s]))
		rand_std_es.append(std_es)

	sub_stes = []
	for s in sub_accs:
		sub_stes.append(stats.sem(sub_accs[s]))

	plot(rand_accs, sub_accs, rand_std_es, sub_stes, args.avg_rands)













