import numpy as np
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt; plt.rcdefaults()
import matplotlib.pyplot as plt
from sklearn import preprocessing
import json
from sklearn.preprocessing import normalize
import sys, argparse

gene_count_dict = dict()
with open("numsamples.json") as f:
    gene_count_dict = json.load(f)

parser = argparse.ArgumentParser(description='Visualization tool for comparing gene expression levels in tissue')

parser.add_argument('--dataset', help='dataset desired to compare from', type=str, required=True)
parser.add_argument('--s1', help='sample class to choose from', type=str, required=True)
parser.add_argument('--s2', help='sample class to choose from', type=str, required=True)
parser.add_argument('--s3', help='sample class to choose from', type=str, required=False)
parser.add_argument('--s4', help='sample class to choose from', type=str, required=False)
parser.add_argument('--stack', help='stacked bar graphs or side-by-side (1 for stacked)', type=int, required=False, default=1)

args = parser.parse_args()

#Example: python compare_means.py "HALLMARK_ADIPOGENESIS" "Bladder" "Lung" will produce BladderLung.png
process = args.dataset
sample1 = args.s1 #"Adipose-Subcutaneous"
sample2 = args.s2#"Adipose-Visceral"
tissues = [sample1,sample2]
maxabsscalar = preprocessing.MaxAbsScaler()

#for Hallmarks, there are 200 genes
average_expressions1 = np.zeros(200)
average_expressions2 = np.zeros(200)
average_expressions = [average_expressions1,average_expressions2]

for i in range(len(tissues)):#go through both samples
    tissue_samples = os.listdir("./datasets/hallmark_subsets/"+process+"/" + tissues[i])#get dat files/samples

    for j in range(gene_count_dict[tissues[i]]):
        #print(tissue_samples[0])
        sample = np.fromfile("./datasets/hallmark_subsets/" + process + "/" + tissues[i]+'/' + tissue_samples[j], dtype=np.float32)
        for d in range(200):
            average_expressions[i][d] += sample[d]

    average_expressions[i] = np.reshape(average_expressions[i], (average_expressions[i].shape[0], 1))

average_expressions = np.hstack((average_expressions[0],average_expressions[1]))
average_expressions = maxabsscalar.fit_transform(average_expressions) 
average_expressions = np.hsplit(average_expressions,2)

if args.stack == 1:
	N = 200
	ind = np.arange(N)    # the x locations for the groups
	width = 0.35       # the width of the bars: can also be len(x) sequence

	p1 = plt.bar(ind, average_expressions[0], width)
	p2 = plt.bar(ind, average_expressions[1], width,
	             bottom=average_expressions[0])

	plt.ylabel('Average Expression Level')
	plt.xlabel('Gene Number')
	plt.title(sample1 + " vs "+ sample2 )
	plt.ylim(0,2)
	plt.legend((p1[0], p2[0]), (sample1, sample2))

	plt.savefig('./graphs/hallmark_subsets/'+ process + '/compares/' + sample1 + sample2 + '.png')
else:
	N = 200
	ind = np.arange(N)    # the x locations for the groups
	width = 0.35       # the width of the bars: can also be len(x) sequence
	opacity = 0.8
	fig, ax = plt.subplots()

	p1 = plt.bar(ind, average_expressions[0], width, alpha=opacity, color='b', label=args.s1)
	p2 = plt.bar(ind + width, average_expressions[1], width, alpha=opacity, color='r', label=args.s2)

	plt.ylabel('Average Expression Level')
	plt.xlabel('Gene Number')
	plt.title(sample1 + " vs "+ sample2 )
	plt.ylim(0,np.max(average_expressions))
	plt.legend((p1[0], p2[0]), (sample1, sample2))
	plt.tight_layout()

	print(args.s1)
	print(str(np.max(average_expressions[0])))
	print(np.argmax(average_expressions[0]))

	print(args.s2)
	print(str(np.max(average_expressions[1])))
	print(np.argmax(average_expressions[1]))

	plt.savefig('./graphs/hallmark_subsets/'+ process + '/compares/' + sample1 + sample2 + '.png')

