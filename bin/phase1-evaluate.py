#!/usr/bin/env python3

'''
This script evaluates the classification potential of gene sets on a dataset.
'''
import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import random
import seaborn as sns
import sklearn.metrics
import sys

import utils



def rotate_xticklabels(angle):
    for tick in plt.gca().get_xticklabels():
        tick.set_horizontalalignment('right')
        tick.set_rotation(angle)



def plot_confusion_matrix(name, y_true, y_pred, classes, output_dir='.'):
    cnf_matrix = sklearn.metrics.confusion_matrix(y_true, y_pred)

    sns.heatmap(cnf_matrix, annot=True, fmt='d', xticklabels=classes, yticklabels=classes)
    plt.tight_layout()
    rotate_xticklabels(45)
    plt.title('Confusion Matrix')
    plt.ylabel('Expected')
    plt.xlabel('Predicted')
    plt.savefig('%s/%s.confusion_matrix.png' % (output_dir, name))
    plt.close()



def evaluate_curated(data, labels, clf, name, genes, n_iters=10, cv=5, n_jobs=None, verbose=True, visualize=False, outfile=None):
    # evaluate gene set n_iters times
    scores = []

    for i in range(n_iters):
        # evaluate gene set
        score, y_true, y_pred = utils.evaluate_gene_set(data, labels, clf, genes, cv=cv, n_jobs=n_jobs)

        scores.append(score)

    # print results
    if verbose:
        n, mu, sigma = len(scores), np.mean(scores), np.std(scores)
        print('%-40s %3d %8.3f %8.3f' % (name, n, mu, sigma))

    # visualize results
    if visualize:
        plot_confusion_matrix(name, y_true, y_pred, classes, output_dir=args.output_dir)

    # write results to output file
    if outfile:
        for score in scores:
            outfile.write('%s\t%.3f\n' % (name, score))



def evaluate_random(data, labels, clf, n_genes, n_iters=100, cv=5, n_jobs=None, verbose=True, outfile=None):
    # evaluate n_iters random sets
    scores = []

    for i in range(n_iters):
        # generate random gene set
        genes = random.sample(list(data.columns), n_genes)

        # evaluate gene set
        score, _, _ = utils.evaluate_gene_set(data, labels, clf, genes, cv=cv, n_jobs=n_jobs)

        scores.append(score)

    # print results
    if verbose:
        n, mu, sigma = len(scores), np.mean(scores), np.std(scores)
        print('%-40s %3d %8.3f %8.3f' % (str(n_genes), n, mu, sigma))

    # write results to output file
    if outfile:
        for score in scores:
            outfile.write('%s\t%.3f\n' % (str(n_genes), score))



if __name__ == '__main__':
    # parse command-line arguments
    parser = argparse.ArgumentParser(description='Evaluate classification potential of gene sets')
    parser.add_argument('--dataset', help='input dataset (samples x genes)', required=True)
    parser.add_argument('--labels', help='list of sample labels', required=True)
    parser.add_argument('--model-config', help='model configuration file (JSON)', required=True)
    parser.add_argument('--model', help='classifier model to use', default='mlp-tf')
    parser.add_argument('--gene-sets', help='list of curated gene sets')
    parser.add_argument('--full', help='Evaluate the set of all genes in the dataset', action='store_true')
    parser.add_argument('--random', help='Evaluate random gene sets', action='store_true')
    parser.add_argument('--random-range', help='range of random gene sizes to evaluate', nargs=3, type=int, metavar=('START', 'STOP', 'STEP'))
    parser.add_argument('--random-iters', help='number of iterations to perform for random classification', type=int, default=100)
    parser.add_argument('--n-jobs', help='number of parallel jobs to use', type=int, default=1)
    parser.add_argument('--cv', help='number of folds for k-fold cross validation', type=int, default=5)
    parser.add_argument('--visualize', help='visualize confusion matrix for each gene set', action='store_true')
    parser.add_argument('--output-dir', help='output directory', default='.')
    parser.add_argument('--outfile', help='output filename', default='phase1-scores.txt')

    args = parser.parse_args()

    # load input data
    print('loading input dataset...')

    df = utils.load_dataframe(args.dataset)
    df_samples = df.index
    df_genes = df.columns

    labels, classes = utils.load_labels(args.labels)

    print('loaded input dataset (%s genes, %s samples)' % (df.shape[1], df.shape[0]))

    # impute missing values
    df.fillna(value=df.min().min(), inplace=True)

    # initialize classifier
    print('initializing classifier...')

    clf = utils.load_classifier(args.model_config, args.model)

    print('initialized %s classifier' % args.model)

    # load gene sets file if it was provided
    if args.gene_sets != None:
        print('loading gene sets...')

        curated_sets = utils.load_gene_sets(args.gene_sets)
        curated_sets = utils.filter_gene_sets(curated_sets, df_genes)

        print('loaded %d gene sets' % (len(curated_sets)))
    else:
        curated_sets = []

    # include the set of all genes if specified
    if args.full:
        curated_sets.append(('FULL', df_genes))

    # initialize list of random set sizes
    if args.random:
        # determine random set sizes from range
        if args.random_range != None:
            print('initializing random set sizes from range...')
            random_sets = range(args.random_range[0], args.random_range[1] + 1, args.random_range[2])

        # determine random set sizes from gene sets
        elif args.gene_sets != None:
            print('initializing random set sizes from curated sets...')
            random_sets = sorted(set([len(genes) for (name, genes) in curated_sets]))

        # print error and exit
        else:
            print('error: --gene-sets or --random-range must be provided to determine random set sizes')
            sys.exit(1)
    else:
        random_sets = []

    print('evaluating gene sets...')

    # initialize output file
    outfile = open(os.path.join(args.output_dir, args.outfile), 'w')
    outfile.write('%s\t%s\n' % ('name', 'score'))

    # evaluate curated gene sets
    for (name, genes) in curated_sets:
        evaluate_curated(df, labels, clf, name, genes, cv=args.cv, n_jobs=args.n_jobs, visualize=args.visualize, outfile=outfile)

    # evaluate random gene sets
    for n_genes in random_sets:
        evaluate_random(df, labels, clf, n_genes, n_iters=args.random_iters, cv=args.cv, n_jobs=args.n_jobs, outfile=outfile)
