#!/usr/bin/env python3

'''
This script uses a random forest classifier to measure the saliency of each
gene in a gene set. Genes with a higher 'importance' are selected as 'candidate'
genes, while the other genes are labeled as 'non-candidate' genes.
'''
import argparse
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import sklearn.ensemble
import sklearn.mixture

import utils



def compute_threshold(genes, scores):
    # fit a Gaussian mixture model to the gene scores
    X = scores.reshape(-1, 1)

    gmm = sklearn.mixture.GaussianMixture(n_components=2)
    gmm.fit(X)

    # compute the intersection between the two modes
    m1 = gmm.means_[0, 0]
    m2 = gmm.means_[1, 0]
    s1 = gmm.covariances_[0, 0, 0]
    s2 = gmm.covariances_[1, 0, 0]

    num = m2*s1**2 - m1*s2**2
    delta = s1 * s2 * math.sqrt((m1 - m2)**2 + 2 * (s1**2 - s2**2) * math.log(s1/s2))
    denom = s1**2 - s2**2

    m = (m1 + m2) / 2
    c1 = (num + delta) / denom
    c2 = (num - delta) / denom

    if abs(c1 - m) < abs(c2 - m):
        threshold = c1
    else:
        threshold = c2

    return threshold



if __name__ == '__main__':
    # parse command-line arguments
    parser = argparse.ArgumentParser(description='Identify candidate / non-candidate genes in a gene set')
    parser.add_argument('--dataset', help='input dataset (samples x genes)', required=True)
    parser.add_argument('--labels', help='list of sample labels', required=True)
    parser.add_argument('--gene-sets', help='list of curated gene sets')
    parser.add_argument('--full', help='Evaluate the set of all genes in the dataset', action='store_true')
    parser.add_argument('--n-jobs', help='number of parallel jobs to use', type=int, default=1)
    parser.add_argument('--threshold', help='manual threshold based on percentile (0-100)', type=float, default=-1)
    parser.add_argument('--visualize', help='visualize candidate threshold', action='store_true')
    parser.add_argument('--output-dir', help='output directory', default='.')

    args = parser.parse_args()

    # load input data
    print('loading input dataset...')

    df = utils.load_dataframe(args.dataset)
    df_samples = df.index
    df_genes = df.columns

    labels, classes = utils.load_labels(args.labels)

    print('loaded input dataset (%s genes, %s samples)' % (df.shape[1], df.shape[0]))

    # load gene sets
    if args.gene_sets != None:
        print('loading gene sets...')

        gene_sets = utils.load_gene_sets(args.gene_sets)
        gene_sets = utils.filter_gene_sets(gene_sets, df_genes)

        print('loaded %d gene sets' % (len(gene_sets)))
    else:
        gene_sets = []

    # include the set of all genes if specified
    if args.full:
        gene_sets.append(('FULL', df_genes))

    # initialize output file
    outfile = open('%s/phase2-rf-genesets.txt' % (args.output_dir), 'w')

    # select candidate genes for each gene set
    for name, genes in gene_sets:
        print('decomposing %s (%d genes)...' % (name, len(genes)))

        # extract dataset
        X = df[genes]
        y = labels

        # normalize dataset
        X = sklearn.preprocessing.MaxAbsScaler().fit_transform(X)

        # determine feature importances
        clf = sklearn.ensemble.RandomForestClassifier(n_estimators=100, n_jobs=args.n_jobs)
        clf.fit(X, y)

        scores = clf.feature_importances_

        # use a percentile threshold if specified
        if args.threshold != -1:
            threshold = np.percentile(scores, args.threshold)

        # otherwise compute threshold automatically
        else:
            threshold = compute_threshold(genes, scores)

        # select candidate genes
        candidate_genes = [gene for i, gene in enumerate(genes) if scores[i] > threshold]

        # plot distribution of gene scores
        if args.visualize:
            sns.distplot(scores)
            ymin, ymax = plt.gca().get_ylim()
            y = [ymin, ymax / 2]
            plt.plot([threshold, threshold], y, 'r')
            plt.title(name)
            plt.savefig('%s/%s-rf-candidate-threshold.png' % (args.output_dir, name))
            plt.close()

        # save results to output file
        outfile.write('\t'.join([name] + candidate_genes) + '\n')
