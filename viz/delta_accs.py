'''

delta_accs.py

	This script allows you to visualize the delta accuracies between random data accuracies
	and a data file of subset accuracies

	The program expects a file that contains the following format:

	Name \t Avg \t StdDev \t Max \t Min
	...
	exname \t 1.0 \t 0.0 \t 1.0 \t1.0
	...

'''



import sys, argparse
import matplotlib.pyplot as plt
import json
import numpy as np

# read in the accuracy files
def read_file(file):
	avg_accs = {}
	with open(file, 'r') as f:
		next(f) # skip the first line which has string header info
		for line in f:
			line = line.split('\t')
			avg_accs[line[0]] = float(line[1])

	return avg_accs


# plot the delta accuracies
def plot(d_accs, gene_counts, out):
	x = np.arange(len(d_accs.keys()))

	difs = [d_accs[s] for s in sorted(d_accs.keys())]

	fig = plt.figure(figsize=(40,20))#figsize=(20,10))


	plt.title("Delta Accuracies") 
	plt.xticks(x, sorted(d_accs.keys()), fontsize=8)
	plt.xticks(rotation=75, ha='right')

	ax1 = fig.add_subplot(111)
	ax2 = fig.add_subplot(111)

	ax1.bar(x, difs, color='#91baff', alpha=0.6)
	ax1.set_ylabel('Delta', color='#91baff')

	ax2 = ax1.twinx()
	ax2.scatter(x, gene_counts, color='#cc0000', s=7)
	ax2.yaxis.tick_right()
	ax2.xaxis.tick_bottom()
	ax2.yaxis.set_label_position('right') 
	ax2.set_ylabel('Gene Count', color='#cc0000')

	plt.gcf().subplots_adjust(bottom=0.3)
	#plt.tight_layout()

	if out:
		fig.savefig(out, bbox_inches='tight')
	else:
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

	# get the delta accs
	gene_counts = []
	delta_accs = {}
	for s in sorted(gene_count_dict.keys()):
		delta = round(sub_accs[s] - rand_accs[str(gene_count_dict[s])], 3)

		if delta < 0:
			delta = 0.0

		delta_accs[s] = delta

		gene_counts.append(gene_count_dict[s])

	plot(delta_accs, gene_counts, args.out)













