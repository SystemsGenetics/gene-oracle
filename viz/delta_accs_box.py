'''

delta_accs_box.py

	This script allows for plotting of the differences between the classification accuracies of random genes vs the
	classification accuracy of hallmark genes

	The program expects the following parameter
	--rand_accs : dir containing random accuracies as a string
	--sub_accs : dir containing subset accuracies as a string
	--sub_count : json with count of genes in subset
	--sub_count_yaxis : optional file to include if you want a different gene list for the y axis than the first subcount
	--out: file to save to


	example: python delta_accs_box.py --rand_accs "../logs/panTCGA_hallmark_random.log" --sub_accs "../logs/panTCGA_hallmark_results.log" --sub_count "../subsets/hallmark_gene_counts.json" --out "test.txt"


'''
import sys, argparse
import matplotlib.pyplot as plt
import json
import numpy as np
import seaborn as sns
import operator as op
import pandas as pd
import collections
import matplotlib.patches as mpatches
import math
from itertools import combinations
from scipy.stats import ttest_ind


def calc_pval(list1,list2):
	#print(list1[0][0])
	for i in range(len(list1)):
		t, p = ttest_ind(list1[i][0], list2[i][0])
		print(p)

def get_mid(list1):
	for i in range(len(list1)):
		print(list1[i][0][2])


# read in the accuracy files
def read_fileAvg(file):
	avg_accs = {}
	with open(file, 'r') as f:
		next(f) # skip the first line which has string header info
		for line in f:
			line = line.split('\t')
			if line[0] not in avg_accs:
				avg_accs[line[0]] = [float(line[1])]
			else:
				avg_accs[line[0]].append(float(line[1]))

	return avg_accs


def read_file(file):
	accs = {}
	with open(file, 'r') as f:
		next(f) # skip the first line which has string header info
		for line in f:
			line = line.split()
			info = []
			info.append(float(line[1]) + float(line[2]))#/math.sqrt(num_samples))# the average plus the standarddev
			info.append(float(line[1]) - float(line[2]))#/math.sqrt(num_samples))# the average minus the standarddev
			info.append(float(line[1]))#avg
			info.append(float(line[3]))#max
			info.append(float(line[4]))#min

			if line[0] not in accs:
				accs[line[0]] = [info]
			else:
				accs[line[0]].append(info)

	return accs


def plotDeltaBoxPlots(ran_accs, sub_accs,out,gene_counts,data_set):

	plt.style.use('ggplot')

	rand_data = ran_accs.values()
	rand_data = list(rand_data)

	sub_data = sub_accs.values()
	sub_data = list(sub_data)

	fig, ax = plt.subplots()
	ax.grid(False)
	fig.set_size_inches(20,10)

	pos = np.arange(len(rand_data))+.13
	ran_bp = ax.boxplot(rand_data,positions=pos,vert =False,widths=.01,whis=0,
                 patch_artist=True, boxprops=dict(color = 'red',edgecolor='red'),medianprops=dict(marker='.',color='k',markersize=4),showfliers=False)
	plt.setp(ran_bp['whiskers'], color="red")

	pos = np.arange(len(sub_data))-.13
	sub_bp = ax.boxplot(sub_data, vert = False,positions=pos,widths=.01,whis=0,
                 patch_artist=True, boxprops=dict(color = 'blue',edgecolor='blue'),medianprops=dict(marker='.',color='k',markersize=4),showfliers=False)
	plt.setp(sub_bp['whiskers'], color="blue")

	random = mpatches.Patch(color='red', label='Random')
	data = mpatches.Patch(color='blue', label=str(data_set))

	ax.legend(handles=[random,data], loc='upper left')

	##Note: This is a custom piece for the panTCGA data and the GTEx data
	if(data_set == "panTCGA"):
		ax.yaxis.set_label_position("right")
		ax.yaxis.tick_right()
		ax.set_yticklabels(gene_counts, fontsize = 8,fontweight='bold', color ='k')
	if(data_set == "GTEx"):
		ax.set_yticklabels(sub_accs.keys(), fontsize = 8,fontweight='bold', color ='k')
	plt.xticks(np.arange(.2, 1, .1))

	plt.tick_params(axis='y')
	ax.set(xlabel="Accuracy", ylabel="Gene Names")
	ax.set_title(str(data_set)+ " vs Random")
	ax.set_aspect(.03)
	plt.tight_layout()
	#plt.show()
	plt.savefig(out)

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Visualize delta accuracy between random genes and subsets')
	parser.add_argument('--rand_accs', help='file random accuracies', type=str, required=True)
	parser.add_argument('--sub_accs', help='file of subset accuracies', type=str, required=True)
	parser.add_argument('--sub_count', help='json with count of genes in subset', type=str, required=True)
	parser.add_argument('--sub_count_yaxis', help='json with count of genes in subset', type=str, required=False)
	parser.add_argument('--out', help='file to save to', type=str, required=False)
	parser.add_argument('--data_set',help='panTCGA or GTEx', type=str, required =True)
	# additional args?
	args = parser.parse_args()

	# read two accuracy files
	rand_accs = read_file(args.rand_accs)
	sub_accs = read_file(args.sub_accs)

	#read a json file containing the count of genes for each subset
	with open(args.sub_count, 'r') as f:
		gene_count_dict = json.load(f)

	if(args.sub_count_yaxis):
		with open(args.sub_count_yaxis, 'r') as f:
			gene_count_dict_y = json.load(f)


	rand_dict = collections.OrderedDict()
	sub_dict = collections.OrderedDict()
	gene_counts = []

	for k,v in sorted(gene_count_dict.items(), key =lambda x: x[1]):
		rand_dict[k] = []
		sub_dict[k] = []
		rand_dict[k].append(rand_accs[str(v)].pop())
		sub_dict[k].append(sub_accs[k][0])
		gene_counts.append(v)
		#print(k)

	#If you specified you wanted another y axis
	if(args.sub_count_yaxis):
		gene_counts = []
		for k,v in sorted(gene_count_dict_y.items(), key =lambda x: x[1]):
			gene_counts.append(v)


	plotDeltaBoxPlots(rand_dict,sub_dict,args.out,gene_counts,args.data_set)
	# for i in gene_counts:
	#  	print(i)
	#calc_pval(rand_dict.values(),sub_dict.values())
	#get_mid(sub_dict.values())
