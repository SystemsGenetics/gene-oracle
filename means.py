import numpy as np
import os
import matplotlib.pyplot as plt; plt.rcdefaults()
import matplotlib.pyplot as plt
from sklearn import preprocessing

tissues = os.listdir("./tissues")
gene_count_dict = dict()

#Create dictionary of the number of samples in each Tissue Directory
for i in range(len(tissues)):#loop through all tissue types
    tissue_samples = os.listdir("./tissues/"+ tissues[i])#get dat files/samples
    for j in range(len(tissue_samples)):
        gene_count_dict[tissues[i]] = len(tissue_samples)

#for Hallmarks, there are 200 genes
average_expressions = np.zeros(200)
for i in range(len(tissues)):#loop through all tissue types
    tissue_samples = os.listdir("./tissues/"+ tissues[i])#get dat files/samples

    for j in range(gene_count_dict[tissues[i]]):
        #print(tissue_samples[0])
        sample = np.fromfile("./tissues/"+ tissues[i]+'/' + tissue_samples[j], dtype=np.float32)
        for d in range(200):
            average_expressions[d] += sample[d]
        
    average_expressions = np.reshape(average_expressions, (average_expressions.shape[0], 1))
    maxabsscalar = preprocessing.MaxAbsScaler()
    average_expressions = maxabsscalar.fit_transform(average_expressions)
    
    #write graph
    plt.clf()#clears plot
    x_range = np.arange(average_expressions.size)
    plt.bar(x_range,average_expressions)
    plt.title(tissues[i])
    plt.ylabel('Expression Level')
    plt.ylim(0,1)
    plt.xticks(x_range)
    plt.show()
    plt.savefig('./images/'+ tissues[i] +'.png')
    
