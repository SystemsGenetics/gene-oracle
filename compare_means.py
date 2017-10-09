import numpy as np
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt; plt.rcdefaults()
import matplotlib.pyplot as plt
from sklearn import preprocessing
import json
from sklearn.preprocessing import normalize

gene_count_dict = dict()
with open("numsamples.json") as f:
    gene_count_dict = json.load(f)

sample1 =  "Adipose-Subcutaneous"
sample2 = "Adipose-Visceral"
tissues = [sample1,sample2]
maxabsscalar = preprocessing.MaxAbsScaler()

#for Hallmarks, there are 200 genes
average_expressions1 = np.zeros(200)
average_expressions2 = np.zeros(200)
average_expressions = [average_expressions1,average_expressions2]

for i in range(len(tissues)):#go through both samples
    tissue_samples = os.listdir("./datasets/hallmark_subsets/HALLMARK_ADIPOGENESIS/"+ tissues[i])#get dat files/samples

    for j in range(gene_count_dict[tissues[i]]):
        #print(tissue_samples[0])
        sample = np.fromfile("./datasets/hallmark_subsets/HALLMARK_ADIPOGENESIS/"+ tissues[i]+'/' + tissue_samples[j], dtype=np.float32)
        for d in range(200):
            average_expressions[i][d] += sample[d]

    average_expressions[i] = np.reshape(average_expressions[i], (average_expressions[i].shape[0], 1))

average_expressions = np.hstack((average_expressions[0],average_expressions[1]))
average_expressions = normalize(average_expressions) 
average_expressions = np.hsplit(average_expressions,2)

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

plt.savefig('./graphs/hallmark_subsets/HALLMARK_ADIPOGENESIS/compares/' + sample1 + sample2 + '.png')
