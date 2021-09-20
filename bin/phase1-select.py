#!/usr/bin/env python3

'''
Use Student's t-test to determine whether a gene set exhibits
significantly higher classification accuracy over equivalent random sets. The
script takes three inputs -- a list of curated sets and their scores, a list of
random sets and their scores, and a list of gene sets and their member
genes -- and produces a list of t-test results for each gene set. The lower the
p-value, the more likely it is that the corresponding gene set performs
significantly better over random.
'''
import argparse
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import operator
import pandas as pd
import scipy.stats

import utils



def plot_delta_boxplots(scores, gene_sets, outfile):
    # gather gene set names and sizes
    names = [name for (name, genes) in gene_sets]
    sizes = [len(genes) for (name, genes) in gene_sets]

    # gather scores into rows for each gene set
    scores_fg = []
    scores_bg = []

    for name, genes in gene_sets:
        name_fg = name
        name_bg = str(len(genes))

        scores_fg.append(scores.loc[scores['name'] == name_fg, 'score'])
        scores_bg.append(scores.loc[scores['name'] == name_bg, 'score'])

    # sort rows by gene set size
    index_array = np.argsort(sizes)
    scores_fg = np.array(scores_fg)[index_array].tolist()
    scores_bg = np.array(scores_bg)[index_array].tolist()

    # create delta box plots
    fig, ax = plt.subplots()

    positions = np.arange(len(scores_fg)) + 0.125
    boxplot_fg = ax.boxplot(
        scores_fg,
        vert=False,
        whis=0,
        positions=positions,
        widths=0.01,
        patch_artist=True,
        showfliers=False,
        boxprops=dict(color='blue', edgecolor='blue'),
        medianprops=dict(marker='.', color='k', markersize=4))
    plt.setp(boxplot_fg['whiskers'], color='blue')

    positions = np.arange(len(scores_bg)) - 0.125
    boxplot_bg = ax.boxplot(
        scores_bg,
        vert=False,
        whis=0,
        positions=positions,
        widths=0.01,
        patch_artist=True,
        showfliers=False,
        boxprops=dict(color='red', edgecolor='red'),
        medianprops=dict(marker='.', color='k', markersize=4))
    plt.setp(boxplot_bg['whiskers'], color='red')

    patch_fg = mpatches.Patch(color='blue', label='Curated')
    patch_bg = mpatches.Patch(color='red', label='Random')

    plt.legend(handles=[patch_fg, patch_bg], loc='upper left')
    plt.xlabel('Accuracy')
    plt.ylabel('Gene Sets')
    plt.xlim(0, 1)
    plt.yticks(np.arange(len(names)), names)
    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()



if __name__ == '__main__':
    # parse command-line arguments
    parser = argparse.ArgumentParser(description='Select gene sets which perform significantly better than equivalent random sets.')
    parser.add_argument('--dataset', help='input dataset (samples x genes)', required=True)
    parser.add_argument('--gene-sets', help='list of curated gene sets', required=True)
    parser.add_argument('--scores', help='list of scores for curated and random gene sets', required=True)
    parser.add_argument('--threshold', help='maximum p-value required for a gene set to be selected', type=float, default=1)
    parser.add_argument('--n-sets', help='maximum number of gene sets that can be selected', type=int, default=-1)
    parser.add_argument('--visualize', help='visualize confusion matrix and ROC for each gene set', action='store_true')
    parser.add_argument('--output-dir', help='output directory', default='.')

    args = parser.parse_args()

    # load input dataset
    df = utils.load_dataframe(args.dataset)
    df_genes = df.columns

    # load gene sets file
    gene_sets = utils.load_gene_sets(args.gene_sets)
    gene_sets = utils.filter_gene_sets(gene_sets, df_genes)

    # load scores file
    scores = pd.read_csv(args.scores, sep='\t')

    # evaluate each curated gene set
    results = []

    for name, genes in gene_sets:
        # perform t-test between gene set scores and background scores
        name_fg = name
        name_bg = str(len(genes))

        t, p = scipy.stats.ttest_ind(
            scores.loc[scores['name'] == name_fg, 'score'],
            scores.loc[scores['name'] == name_bg, 'score'],
            equal_var=False)

        results.append((name, genes, p))

        # print result
        print('%-40s %0.12f' % (name, p))

    # visualize gene sets if specified
    if args.visualize:
        plot_delta_boxplots(scores, gene_sets, '%s/phase1-genesets.png' % (args.output_dir))

    # select gene sets
    results.sort(key=operator.itemgetter(2))
    results = [(name, genes) for (name, genes, p) in results if p <= args.threshold]

    if args.n_sets != -1:
        results = results[0:args.n_sets]

    # save selected gene sets to output file
    outfile = open('%s/phase1-genesets.txt' % (args.output_dir), 'w')

    for name, genes in results:
        outfile.write('\t'.join([name] + genes) + '\n')
