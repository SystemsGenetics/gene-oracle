'''
#Todo: compare light vs dark

heatmap.py

	This script allows you to create heatmaps/clustermaps the frequency of genes in combos that persist in the
    pool of combos that are run through the network

	The program expects a directory path that contains the following formated files:

    graph_type : clustermap or heatmap
	analysis: all, topten
    directory: /...path.../
    num_genes: #
    hallmark:
    dataset: gtex or panTCGA

    halmark file format: 'hallmark_#_gene_accuracy.txt'
                        Contents: [list of genes] accuracy


example:  python heatmap_gen.py --graph_type "heatmap" --analysis all "--directory" "../logs/hallmark_hedgehog_signaling_panTCGA/" --num_genes 36 --hallmark "hallmark_hedgehog_signaling" --dataset "panTCGA"

'''

import sys, argparse
import matplotlib.pyplot as plt
import json
import numpy as np
import seaborn as sns
import pandas as pd


# takes the data in the log files for the specific hallmark subset
# and loads it into List of Dataframes with columns of "Genes" and "Accuracy"
def getDataFromLog(directory,num_genes,hallmark,dataset):
	comboList = []
	#process datafiles
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

	return comboList

#Load panTCGA Data
def loadpanTCGA(directory,genes):
    ensembles = np.loadtxt(directory + "gene_list.txt",dtype=str)

	#clean
    for i in range(len(ensembles)):
        ensembles[i] = ensembles[i].replace("b'","")
        ensembles[i] = ensembles[i].replace("'","")

	#set initial genes in row one
    for gene in ensembles:
        genes.loc[1,gene] = 1

    return genes

#Load GTEx Data
def loadGTEx(hallmark,genes):
    with open('../subsets/gene_dict.json', 'r') as fp:
            gene_counts = json.load(fp)

	#clean and set initial genes in row one
    for gene in gene_counts[hallmark.upper()]:
         gene = gene.replace("('","")
         gene = gene.replace("')","")
         gene =  gene.replace("'","")
         genes.loc[1,gene] = 1

    return genes

# Generates a dataframe for the dataset specified
def prepareDataSetAll(directory,num_genes,hallmark,dataset):

    genes = pd.DataFrame(index =range(1,num_genes))

	#Choose Dataset
    if(dataset == 'GTEx'):#GTEx
        genes = loadGTEx(hallmark,genes)

	#if panTCGA is chosen, then the genes must be converted to ensembles
    if(dataset == "panTCGA"):
        genes = loadpanTCGA(directory,genes)

	#normalize the first row
    genes.loc[1] = genes.loc[1].div(num_genes)
    return genes

# Fills out the genes dataframe with the frequency of a gene in all the combos
def freqCountAll(directory,num_genes,hallmark,dataset):
    comboList = getDataFromLog(directory,num_genes,hallmark,dataset)
    genes = prepareDataSetAll(directory,num_genes,hallmark,dataset)

	#Count the frequency of genes in the combos
    for i in range(1, num_genes):
        print("Counting"+ str(i))
        count = 0
        for combo in comboList[i-1]['Genes']:
             if(i == 1):#special case for the first gene, need to fix but works
                array = []
                array.append(eval(combo))
             else:#This is for the rest of the genes
                array = eval(combo)
                array = list(array)
			#Count Genes
             for gene in genes.columns.values:
                 if gene in array:
                     if np.isnan(genes.loc[i,gene]):
                         genes.loc[i,gene] = 1
                     else:
                         genes.loc[i,gene] += 1
        count = genes.loc[i].sum()
        genes.loc[i] = genes.loc[i].div(count)

    if(dataset == "panTCGA"):#need to convert for readability
        with open('../data/ensembles_to_hallmark_id.json', 'r') as fp:
            ensembles_list = json.load(fp)
        for i in range(num_genes):
            genes = genes.rename(columns={genes.columns[i]:ensembles_list[genes.columns[i]]})

    return genes

def prepareDataSetTopTen(directory,num_genes,hallmark,dataset):
    genes = pd.DataFrame(index =range(1,num_genes))

	#Choose Dataset
    if(dataset == 'GTEx'):#GTEx
        genes = loadGTEx(hallmark,genes)

	#if panTCGA is chosen, then the genes must be converted to ensembles
    if(dataset == "panTCGA"):
        genes = loadpanTCGA(directory,genes)

    return genes

def freqCountTopTen(directory,num_genes,hallmark,dataset):
    comboList = getDataFromLog(directory,num_genes,hallmark,dataset)
    genes = prepareDataSetTopTen(directory,num_genes,hallmark,dataset)

	#Count the frequency of genes in the combos
    for i in range(1, num_genes):
        print("Counting"+ str(i))
        count = 0
        best_combo_list = comboList[i-1].sort_values("Accuracy",ascending=False)[0:10]
        genes.loc[i] = 0
        print(best_combo_list)
        for combo in best_combo_list["Genes"]:
             if(i == 1):
                gene_array = []
                gene_array.append(eval(combo))
             else:
                gene_array = eval(combo)
                gene_array = list(gene_array)

             for gene in gene_array:
                if genes.loc[i,str(gene)] == 0:
                   genes.loc[i,str(gene)] = 1
                else:
                   genes.loc[i,str(gene)] += 1
        count = genes.loc[i].sum()
        genes.loc[i] = genes.loc[i].div(count)

    if(dataset == "panTCGA"):#need to convert for readability
        with open('../data/ensembles_to_hallmark_id.json', 'r') as fp:
            ensembles_list = json.load(fp)
        for i in range(num_genes):
            genes = genes.rename(columns={genes.columns[i]:ensembles_list[genes.columns[i]]})
    return genes

#Generates a plot
def plot(graph_type,hallmark,data,analysis,dataset):
    print("Plotting..")
    sns.reset_defaults()
    fig, ax = plt.subplots(figsize=(20,10))
    if(graph_type == 'clustermap'):
        sns.clustermap(data=data,cmap="Blues",linewidths=1, linecolor='black')
    if(graph_type == 'heatmap'):
        sns.heatmap(data=data,cmap="Blues",linewidths=1, linecolor='black',vmin=0, vmax=.2)

    plt.xticks(rotation=90)
    plt.yticks(rotation=0)
    plt.xlabel("Genes")
    plt.ylabel("Combo")
    plt.title(hallmark +"_"+ dataset +"_"+ analysis)
    plt.savefig("../graphs/"+hallmark + dataset+ str(analysis) + '.png')



if __name__ == '__main__':

        parser = argparse.ArgumentParser(description='Visualize Gene Frequency')
        parser.add_argument('--graph_type', help='type of graph', type=str, required=True)
        parser.add_argument('--directory', help='location of files', type=str, required=True)
        parser.add_argument('--analysis', help='type of data analysis', type=str, required=True)
        parser.add_argument('--num_genes',help="number of genes in hallmark", type=int, required=True)
        parser.add_argument('--hallmark',help="hallmark subset", type=str, required=True)
        parser.add_argument('--dataset',help="data", type=str, required=True)


        args = parser.parse_args()
        if(str(args.analysis) == "all"):
             data = freqCountAll(args.directory,args.num_genes,args.hallmark,args.dataset)
        if(str(args.analysis) == "topten"):
             data = freqCountTopTen(args.directory,args.num_genes,args.hallmark,args.dataset)
        print(data)
        #plot(args.graph_type,args.hallmark,data,args.analysis,args.dataset)
