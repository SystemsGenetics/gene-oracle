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
import numpy as np
import pandas as pd


def read_file(subset,dataset,set):
    with open("../lists.json") as f:
        data = json.load(f)

    return data[subset][dataset][set]

def count_interactions(gene_list):
    total_interactions = 0
    #read in BioGRID
    biogrid = pd.read_table("../data/BIOGRID-ALL-3.5.165.tab2.txt", sep='\t')

    # only use human genes
    biogrid = biogrid[biogrid['Organism Interactor A'] == 9606]
    biogrid = biogrid[biogrid['Organism Interactor B'] == 9606]

    pairs = biogrid[['Official Symbol Interactor A','Official Symbol Interactor B']].values

    for gene in gene_list:
        gene_pairs = []

        # get locations where gene occurs in pairs
        locs = np.where(pairs==gene)

        # get list of genes that occur in interactions... remove gene of interest
        uniq = list(np.unique(pairs[locs[0]]))
        if(len(uniq)):
            uniq.remove(str(gene))

        pairs = np.delete(pairs, locs[0], axis=0)

        total_interactions += len(uniq)

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

    if(args.set == "RANDOM"):
        print("Running Randoms...")
        total_ran_interactions = 0
        for i in range(5):
            total_ran_interactions += count_interactions(gene_list[i])
        print(total_ran_interactions / 5)#average

    else:

        print(count_interactions(gene_list))
