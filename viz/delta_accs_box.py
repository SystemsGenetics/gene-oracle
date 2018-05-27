'''

delta_accs_box.py

	This script allows for plotting of the differences between the classification accuracies of random genes ve the
	classification accuracy of hallmark genes

	The program expects the following parameter
	--rand_accs : dir containing random accuracies as a string
	--sub_accs : dir containing subset accuracies as a string
	--sub_count : json with count of genes in subset
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


def read_file(file,num_samples):
	accs = {}
	with open(file, 'r') as f:
		next(f) # skip the first line which has string header info
		for line in f:
			line = line.split()
			info = []
			info.append(float(line[1]) + float(line[2])/math.sqrt(num_samples))# the average plus the standarddev
			info.append(float(line[1]) - float(line[2])/math.sqrt(num_samples))# the average minus the standarddev
			info.append(float(line[1]))#avg
			info.append(float(line[3]))#max
			info.append(float(line[4]))#min

			if line[0] not in accs:
				accs[line[0]] = [info]
			else:
				accs[line[0]].append(info)

	return accs


# plot the delta accuracies
def plotDelta(d_accs, gene_counts, out):
	rc = {'font.size' : 6}
	sns.set(rc)

	data = d_accs.values()
	fixed_data = list(data)

	fig, ax = plt.subplots()
	fig.set_size_inches(20,10)

	ax = sns.boxplot(data=fixed_data,orient='h')
	ax.set_yticklabels(delta_accs.keys())
	ax.set(xlabel="Delta Accuracy", ylabel="Gene Name")
	ax.set_title("panTCGA Delta Accuracy")


	plt.show()

def plotDouble(ran_accs, sub_accs,out,gene_counts,data_set):

	plt.style.use('ggplot')


	rand_data = ran_accs.values()
	rand_data = list(rand_data)

	sub_data = sub_accs.values()
	sub_data = list(sub_data)

	print(rand_data)
	print(sub_data)
	#print(str(len(sub_data)) +" " +str(type(sub_data)))
	#print(str(len(rand_data)) + " " + str(type(rand_data)))

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

	if(data_set == "panTCGA"):
		ax.yaxis.set_label_position("right")
		ax.yaxis.tick_right()
		ax.set_yticklabels(gene_counts, fontsize = 8,fontweight='bold', color ='k')
	else:
		ax.set_yticklabels(sub_accs.keys(), fontsize = 8,fontweight='bold', color ='k')
	plt.xticks(np.arange(.2, 1, .1))

	plt.tick_params(axis='y')
	ax.set(xlabel="Accuracy", ylabel="Gene Names")
	ax.set_title(str(data_set)+ " vs Random")
	ax.set_aspect(.03)
	plt.tight_layout()

	plt.savefig(out)



if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Visualize delta accuracy between random genes and subsets')
	parser.add_argument('--rand_accs', help='file random accuracies', type=str, required=True)
	parser.add_argument('--sub_accs', help='file of subset accuracies', type=str, required=True)
	parser.add_argument('--sub_count', help='json with count of genes in subset', type=str, required=True)
	parser.add_argument('--out', help='file to save to', type=str, required=False)
	parser.add_argument('--gene_count', help='gene counts file', type=str, required=False)
	parser.add_argument('--data_set',help='panTCGA or GTEx', type=str, required =True)
	# additional args?
	args = parser.parse_args()

	# read two accuracy files
	rand_accs = read_file(args.rand_accs,50)
	sub_accs = read_file(args.sub_accs,10)

	# print (rand_accs)

	# read a json file containing the count of genes for each subset
	with open(args.sub_count, 'r') as f:
		gene_count_dict = json.load(f)


	gene_counts = []
	rand_dict = collections.OrderedDict()
	sub_dict = collections.OrderedDict()
	gene_counts = []

	for k,v in sorted(gene_count_dict.items(), key =lambda x: x[1]):

		rand_dict[k] = []
		sub_dict[k] = []


		rand_dict[k].append(rand_accs[str(v)].pop())
		sub_dict[k].append(sub_accs[k][0])#append(round(sub_accs[k][i] - current_gene_count[i], 3))

		#print(str(gene_count_dict[k]) +" : " + str(k) + " : " + str(rand_dict[k]))
		#print(str(gene_count_dict[k]) +" : " + str(k) + " : " + str(sub_dict[k]))

	gene_counts = []
	with open(args.gene_count, 'r') as fp:
		for line in fp:
			line = line.split("\n")
			gene_counts.append(line[0])


	print(gene_counts)



	# sanity check
	# for key, value in delta_accs.items():
	# 	print(str(key) + " , " + str(value) )



	plotDouble(rand_dict,sub_dict,args.out,gene_counts,args.data_set)
