'''
    This script can be used to find the number unique interactions in a list of genes based on the BIOGRID database.

    To get the data use: wget https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-3.5.165/BIOGRID-ALL-3.5.165.tab2.zip

    Arguments:
        - subset in the config file: hedgehog, angiogenesis, notch, mycv2
        - dataset in the config file: GTEx or TCGA
        - set in the config file: CAN, NONCAN, RANDOM

'''

import sys, argparse
import json
import csv
import pandas as pd


def read_file(subset,dataset,set):
    with open("lists.json") as f:
        data = json.load(f)

    return data[subset][dataset][set]

def count_interactions(gene_list):
    gene_dict = {}
    total_interactions = 0
    #read in BioGRID
    biogrid = pd.read_table("./BIOGRID-ALL-3.5.165.tab2.txt", sep='\t',usecols=[7,8])
    pairs = biogrid[['Official Symbol Interactor A','Official Symbol Interactor B']].values

    for gene in gene_list:
        #print("Gene: " + str(gene))
        current_iterations = 0
        gene_pairs = []
        #get all pairs
        for pair in pairs:
            if gene in pair:
                gene_pairs.append(pair)
        #get list of all pairs
        pairs_list = []
        [pairs_list.append(g) for pair in gene_pairs for g in pair]
        uniq = set(pairs_list)
        uniq.remove(str(gene))
        gene_dict[str(gene)] = uniq

        current_iterations = len(uniq)#-1 for the gene
        #check if prev interaction recorded
        for new_gene in uniq:
            if new_gene in gene_dict.keys():
                for old_gene in gene_dict[new_gene]:
                    if old_gene == gene:
                        print("repeat")
                        current_iterations += -1
        #add genes to the total_interactions ( minus 1 for the gene searched)
        #print(str(gene) + "interactions is " + str(current_iterations))
        total_interactions += current_iterations
        #print("total_interactions is " + str(total_interactions))
    return total_interactions

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Get the Interactions of the given gene set')
    parser.add_argument('--subset', help='Hallmark Subset', type=str, required=True)
    parser.add_argument('--dataset', help='DataSet Type', type=str, required=True)
    parser.add_argument('--set', help='Set Within the Subset', type=str, required=True)

    args = parser.parse_args()
    gene_list = read_file(args.subset,args.dataset,args.set)
    print("Current Gene List is ...")
    print(gene_list)

    if(args.set == "RANDOM"){

    }

    print(count_interactions(gene_list))
