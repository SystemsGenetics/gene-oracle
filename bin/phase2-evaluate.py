#!/usr/bin/env python3

'''
Decompose a gene set into subsets and evaluates them in order to
identify the subsets with the highest classification potential.
'''
import argparse
import itertools
import numpy as np
import operator
import os
import pandas as pd
import random

import utils



def load_scores(filename):
    infile = open(filename, 'r')
    lines = [line.strip() for line in infile]
    lines = [line.split('\t') for line in lines]
    subsets = [(line[0].split(','), float(line[1])) for line in lines]

    return subsets



def save_scores(filename, subsets):
    outfile = open(filename, 'w')

    for subset, score in subsets:
        outfile.write('%s\t%0.3f\n' % (','.join(subset), score))



def select_subsets(prev_subsets, genes, n_subsets=50, r=0.5):
    # sort previous subsets by score (descending)
    prev_subsets.sort(key=operator.itemgetter(1), reverse=True)
    prev_subsets = [s[0] for s in prev_subsets]

    # select the highest scoring subsets from prev subsets
    seed_subsets = prev_subsets[0:n_subsets]
    prev_subsets = prev_subsets[n_subsets:]

    # additionally select random subsets from the remaining prev subsets
    n_random = min(int(r * n_subsets), len(prev_subsets))
    seed_subsets += random.sample(prev_subsets, n_random)

    # generate new subsets by augmenting the seed subsets with individual genes
    subsets = []

    for seed_subset in seed_subsets:
        # determine the set of genes not in the seed subset
        extra_genes = list(set(genes) - set(seed_subset))

        # generate new subsets by appending each extra gene to seed subset
        subsets += [(seed_subset + [gene]) for gene in extra_genes]

    # remove duplicate sets
    subsets = [sorted(subset) for subset in subsets]
    subsets = [list(s) for s in set(tuple(s) for s in subsets)]

    return subsets



def chunk_select(genes, k, infile=None):
    n_genes = len(genes)

    # generate all combinations of size k
    if k <= 3 or n_genes - k <= 1:
        subsets = [list(s) for s in itertools.combinations(genes, k)]

    # or select some combinations using a heuristic
    elif infile != None:
        # load subsets from previous iteration
        prev_subsets = load_scores(infile)

        # generate new subsets of size k from previous subsets
        subsets = select_subsets(prev_subsets, genes)

    # augment with empty scores
    subsets = [(subset, 0) for subset in subsets]

    return subsets



def chunk_evaluate(df, labels, clf, subsets, outfile):
    # evaluate each subset
    results = [utils.evaluate_gene_set(df, labels, clf, genes, n_jobs=args.n_jobs) for genes, _ in subsets]
    subsets = [(genes, score) for ((genes, _), (score, _, _)) in zip(subsets, results)]

    # save results to output file
    save_scores(outfile, subsets)



if __name__ == '__main__':
    # parse command-line arguments
    parser = argparse.ArgumentParser(description='Generate and evaluate subsets of a gene set.')
    parser.add_argument('--dataset', help='input dataset (samples x genes)', required=True)
    parser.add_argument('--labels', help='list of sample labels', required=True)
    parser.add_argument('--gene-sets', help='list of curated gene sets')
    parser.add_argument('--model-config', help='model configuration file (JSON)', required=True)
    parser.add_argument('--model', help='classifier model to use', default='mlp-tf')
    parser.add_argument('--random', help='Evaluate random gene sets', action='store_true')
    parser.add_argument('--random-range', help='range of random gene sizes to evaluate', nargs=2, type=int)
    parser.add_argument('--n-jobs', help='number of parallel jobs to use', type=int, default=1)
    parser.add_argument('--logdir', help='directory where logs are stored', required=True)
    parser.add_argument('--chunk-geneset', help='current gene set for chunk runs')
    parser.add_argument('--chunk-iteration', help='current iteration for chunk runs', type=int)
    parser.add_argument('--chunk-op', help='operation to perform for chunk runs', choices=['select', 'evaluate'])
    parser.add_argument('--chunk-infile', help='input file for chunk runs')
    parser.add_argument('--chunk-outfile', help='output file for chunk runs')

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

        gene_sets = utils.load_gene_sets(args.gene_sets)
        gene_sets = utils.filter_gene_sets(gene_sets, df_genes)

        print('loaded %d gene sets' % (len(gene_sets)))
    else:
        gene_sets = []

    # generate random gene sets if specified
    if args.random:
        # determine random set sizes from range
        if args.random_range != None:
            print('initializing random set sizes from range...')
            random_sets = range(args.random_range[0], args.random_range[1] + 1)

        # determine random set sizes from gene sets
        elif args.gene_sets != None:
            print('initializing random set sizes from gene sets...')
            random_sets = sorted(set([len(genes) for (name, genes) in gene_sets]))

        # print error and exit
        else:
            print('error: --gene-sets or --random-range must be provided to determine random set sizes')
            sys.exit(1)

        # generate random gene sets
        for n_genes in random_sets:
            name = 'random-%d' % n_genes
            genes = random.sample(list(df_genes), n_genes)

            gene_sets.append((name, genes))

    # initialize log directory
    os.makedirs(args.logdir, exist_ok=True)

    # perform chunk run if specified
    if args.chunk_geneset != None and args.chunk_iteration != None:
        print()
        print('performing iteration %d for gene set %s' % (args.chunk_iteration, args.chunk_geneset))

        # search gene set from list
        genes = next(genes for (name, genes) in gene_sets if name == args.chunk_geneset)

        # perform selection or evaluation for the current iteration
        if args.chunk_op == 'select':
            subsets = chunk_select(genes, args.chunk_iteration, args.chunk_infile)

            save_scores(args.chunk_outfile, subsets)

        elif args.chunk_op == 'evaluate':
            subsets = load_scores(args.chunk_infile)

            chunk_evaluate(df, labels, clf, subsets, args.chunk_outfile)

    # otherwise perform full combinatorial analysis on each gene set
    else:
        for name, genes in gene_sets:
            print()
            print('decomposing %s (%d genes)...' % (name, len(genes)))

            # perform combinatorial analysis
            for k in range(1, len(genes) + 1):
                # select subsets
                print('  selecting subsets of size %d' % k)

                subsets = chunk_select(genes, k, '%s/%s_scores_%03d.txt' % (args.logdir, name, k - 1))

                # evaluate subsets
                print('  evaluating %d subsets...' % len(subsets))

                chunk_evaluate(df, labels, clf, subsets, '%s/%s_scores_%03d.txt' % (args.logdir, name, k))
