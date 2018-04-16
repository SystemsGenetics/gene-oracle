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
import collections

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

def plotDouble(ran_accs, sub_accs,out):
	
	rc = {'font.size' : 6}
	sns.set(rc)


	rand_data = ran_accs.values()
	#fixed_data_rand = list(rand_data)

	sub_data = sub_accs.values()
	#fixed_data_sub = list(sub_data)

	print(len(sub_data))
	print(len(rand_data))

	fig, ax = plt.subplots()
	fig.set_size_inches(20,10)

	pos = np.arange(len(rand_data))+.2 
	ran_bp = ax.boxplot(rand_data,positions=pos,vert =False,widths=.25, 
                 patch_artist=True, boxprops=dict(facecolor="red",edgecolor='red',alpha = .5))
	plt.setp(ran_bp['whiskers'], color="red")

	pos = np.arange(len(sub_data))-.2 
	sub_bp = ax.boxplot(sub_data, vert = False,positions=pos,widths=.25, 
                 patch_artist=True, boxprops=dict(edgecolor='blue',alpha=.3))
	plt.setp(sub_bp['whiskers'], color="blue")

	ax.legend([ran_bp["boxes"][0], sub_bp["boxes"][0]], ['random', 'hallmark'], loc='upper right')


	#ax = sns.boxplot(data=fixed_data_rand,orient='h',  color= 'blue')
	#ax = sns.boxplot(data=fixed_data_sub,orient='h', color ='red')
	ax.set_yticklabels(sub_accs.keys())
	plt.tick_params(axis='y')
	ax.set(xlabel="Accuracy", ylabel="Gene Name")
	ax.set_title("Hallmark vs Random")
	ax.set_aspect(.009)
	plt.tight_layout()

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
	sub_accs = read_file(args.sub_accs)

	# read a json file containing the count of genes for each subset
	with open(args.sub_count, 'r') as f:
		gene_count_dict = json.load(f)


	#sort based on gene count

	#print(gene_count_dict)


	# get the delta accs
	gene_counts = []
	delta_accs = collections.OrderedDict()
	for k,v in sorted(gene_count_dict.items(), key =lambda x: x[1]):
		delta_accs[k] = []
		current_gene_count = rand_accs[str(v)]

		for i in range(len(current_gene_count)):
			delta_accs[k].append(round(sub_accs[k][i] - current_gene_count[i], 3))


		print(str(gene_count_dict[k]) +" : " + str(k) + " : " + str(delta_accs[k]))

	##plotDelta(delta_accs, gene_counts, args.out)

	gene_counts = []
	rand_dict = collections.OrderedDict()
	sub_dict = collections.OrderedDict()
	for k,v in sorted(gene_count_dict.items(), key =lambda x: x[1]):
		rand_dict[k] = []
		sub_dict[k] = []
		current_gene_count = rand_accs[str(v)]

		for i in range(len(current_gene_count)):
			rand_dict[k].append(current_gene_count[i])
			sub_dict[k].append(sub_accs[k][i])#append(round(sub_accs[k][i] - current_gene_count[i], 3))


		print(str(gene_count_dict[k]) +" : " + str(k) + " : " + str(rand_dict[k]))	
		print(str(gene_count_dict[k]) +" : " + str(k) + " : " + str(sub_dict[k]))	




	# sanity check
	# for key, value in delta_accs.items():
	# 	print(str(key) + " , " + str(value) )


	
	plotDouble(rand_dict,sub_dict,args.out)