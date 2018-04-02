'''

heatmap.py

	This script allows you to create heatmaps/clustermaps the frequency of genes in combos that persist in the
    pool of combos that are run through the network

	The program expects a directory path that contains the following formated files:

    graph_type : clustermap or heatmap
    directory: /...path.../
    num_genes: #
    hallmark:
    dataset: gtex or panTCGA

    halmark file format: 'hallmark_#_gene_accuracy.txt'
                        Contents: [list of genes] accuracy


example:  python heatmap_gen.py --graph_type "heatmap" "--directory" "../logs/hallmark_hedgehog_signaling_panTCGA/" --num_genes 36 --hallmark "hallmark_hedgehog_signaling" --dataset "panTCGA"

'''



import sys, argparse
import matplotlib.pyplot as plt
import json
import numpy as np
import seaborn as sns
import pandas as pd

genes = pd.DataFrame()
comboList = []

def process(directory,num_genes,hallmark,dataset):
    global genes
    global comboList
    for i in range(1,num_genes):
        print("Processing..." + str(i))
        combos = pd.DataFrame( columns = ["Genes", "Accuracy"])
        count = 0
        with open(directory + hallmark+"_" + str(i)+'_gene_accuracy.txt' , "r") as f:
            for line in f:
                (key, val) = line.split("\t")
                val =val.replace("\n", "")
                key = key.replace(",)",")")
                combos.loc[count,"Genes"] = key
                combos.loc[count,"Accuracy"] = val
                count += 1
            comboList.append(combos)

    genes = pd.DataFrame(index =range(1,num_genes))

    if(dataset == 'gtex'):#Gtex
        with open('../subsets/gene_dict.json', 'r') as fp:
                gene_counts = json.load(fp)

        for gene in gene_counts[hallmark.upper()]:
             gene = gene.replace("('","")
             gene = gene.replace("')","")
             gene =  gene.replace("'","")
             genes.loc[1,gene] = 1
             #genes.loc[2,gene] = 1
    else:#panTCGA
        ensembles = np.loadtxt(directory + "gene_list.txt",dtype=str)
        with open('../data/ensembles_to_hallmark_id.json', 'r') as fp:
                ensembles_list = json.load(fp)

        for i in range(len(ensembles)):
            ensembles[i] = ensembles[i].replace("b'","")
            ensembles[i] = ensembles[i].replace("'","")

        for gene in ensembles:
             genes.loc[1,ensembles_list[gene]] = 1
             #genes.loc[2,ensembles_list[gene]] = 1


    genes.loc[1] = genes.loc[1].div(num_genes)
    #genes.loc[2] = genes.loc[2].div(num_genes)

    for i in range(1, num_genes):
        print("Counting"+ str(i))
        count = 0
        for combo in comboList[i-1]['Genes']:
             if(i == 1):
                array = []
                array.append(eval(combo))
             else:
                array = eval(combo)
                array = list(array)
             #break
             #array = list(array)
             if(dataset != 'gtex'):
                 for j in range(len(array)):
                     array[j] = ensembles_list[array[j]]
             #need to convert here
             for gene in genes.columns.values:
                 if gene in array:
                     if np.isnan(genes.loc[i,gene]):
                         genes.loc[i,gene] = 1
                     else:
                         genes.loc[i,gene] += 1
        count = genes.loc[i].sum()
        genes.loc[i] = genes.loc[i].div(count)


def plot(graph_type,hallmark,dataset):
    print("Plotting..")
    sns.reset_defaults()
    fig, ax = plt.subplots(figsize=(20,10))
    if(graph_type == 'clustermap'):
        sns.clustermap(data=genes,cmap='viridis',linewidths=1, linecolor='yellow')
    if(graph_type == 'heatmap'):
        sns.heatmap(data=genes,cmap='viridis',linewidths=1, linecolor='yellow')

    plt.xticks(rotation=90)
    plt.yticks(rotation=0)
    plt.xlabel("Genes")
    plt.ylabel("Combo")
    plt.title(hallmark + dataset)
    plt.savefig("../graphs/"+hallmark + dataset+ '.png')



if __name__ == '__main__':

        parser = argparse.ArgumentParser(description='Visualize Gene Frequency')
        parser.add_argument('--graph_type', help='type of graph', type=str, required=True)
        parser.add_argument('--directory', help='location of files', type=str, required=True)
        parser.add_argument('--num_genes',help="number of genes in hallmark", type=int, required=True)
        parser.add_argument('--hallmark',help="hallmark subset", type=str, required=True)
        parser.add_argument('--dataset',help="data", type=str, required=True)


        args = parser.parse_args()
        #graph_type = read_file(args.graph_type)
        #directory  = read_file(args.directory)
        #num_genes  = read_file(args.num_genes)
        #hallmark   = read_file(args.hallmark)

        process(args.directory, args.num_genes, args.hallmark,args.dataset)
        plot(args.graph_type, args.hallmark, args.dataset)
#flip the accuracy polt
#boxplots
