#!/usr/bin/env python3

import argparse
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import random
import sklearn.datasets
import sklearn.manifold
import sys

import utils



if __name__ == '__main__':
    # parse command-line arguments
    parser = argparse.ArgumentParser(description='Create a synthetic classification dataset')
    parser.add_argument('--n-samples', help='number of samples', type=int, default=100)
    parser.add_argument('--n-genes', help='number of genes', type=int, default=20)
    parser.add_argument('--n-classes', help='number of classes', type=int, default=2)
    parser.add_argument('--n-sets', help='number of gene sets', type=int, default=10)
    parser.add_argument('--dataset', help='name of dataset file', default='example.emx.txt')
    parser.add_argument('--labels', help='name of label file', default='example.labels.txt')
    parser.add_argument('--gene-sets', help='name of gene sets file', default='example.genesets.txt')
    parser.add_argument('--visualize', help='create t-SNE plot of dataset', action='store_true')

    args = parser.parse_args()

    # create synthetic dataset
    n_informative = args.n_genes // 10
    n_redundant = args.n_genes - n_informative

    x, y = sklearn.datasets.make_classification(
        args.n_samples,
        args.n_genes,
        n_informative=n_informative,
        n_redundant=n_redundant,
        n_classes=args.n_classes,
        n_clusters_per_class=1)

    # initialize class names
    classes = ['class-%02d' % i for i in range(args.n_classes)]
    y = [classes[y_i] for y_i in y]

    # initialize gene names, sample names
    x_samples = ['sample-%08d' % i for i in range(args.n_samples)]
    x_genes = ['gene-%06d' % i for i in range(args.n_genes)]

    # initialize dataframes
    x = pd.DataFrame(x, index=x_samples, columns=x_genes)
    y = pd.DataFrame(y, index=x_samples)

    # create synthetic gene sets
    gene_sets = []

    for i in range(args.n_sets):
        n_genes = random.randint(5, min(max(10, args.n_genes // 10), args.n_genes))
        genes = random.sample(x_genes, n_genes)

        gene_sets.append(['gene-set-%03d' % i] + genes)

    # visualize dataset if specified
    if args.visualize:
        # compute t-SNE embedding
        x_tsne = sklearn.manifold.TSNE().fit_transform(x)

        # plot t-SNE embedding by class
        fig, ax = plt.subplots()
        colors = cm.rainbow(np.linspace(0, 1, len(classes)))

        for c in classes:
            indices = (y[0] == c)
            ax.scatter(x_tsne[indices, 0], x_tsne[indices, 1], label=c, alpha=0.75)

        plt.subplots_adjust(right=0.75)
        ax.set_axis_off()
        ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
        plt.savefig('%s.tsne.png' % (args.dataset.split('.')[0]))
        plt.close()

    # save dataset to file
    utils.save_dataframe(args.dataset, x)

    # save labels to file
    y.to_csv(args.labels, sep='\t', header=None)

    # save gene sets to file
    f = open(args.gene_sets, 'w')
    f.write('\n'.join(['\t'.join(gene_set) for gene_set in gene_sets]))
