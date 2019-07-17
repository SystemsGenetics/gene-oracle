#!/bin/bash
# Example usage of gene-oracle on a synthetic dataset.

# create synthetic input data
python bin/make-input-data.py \
	--n-samples 1000 \
	--n-genes 20 \
	--n-classes 2 \
	--n-sets 20

# evaluate gene sets
python bin/phase1-evaluate.py \
	--dataset      example_data.txt \
	--labels       example_labels.txt \
	--model-config example/models.json \
	--model        lr \
	--gene-sets    example_genesets.txt \
	--random \
	--random-iters 10 \
	--cv           5 \
	--n-jobs       1 \
	--outfile      phase1-scores.txt

# select gene sets which score higher over random
python bin/phase1-select.py \
	--scores    phase1-scores.txt \
	--gene-sets example_genesets.txt \
	--threshold 1 \
	--n-sets 5 \
	--outfile phase1-genesets.txt

# perform combinatorial analysis on selected gene sets
python bin/phase2-evaluate.py \
	--dataset      example_data.txt \
	--labels       example_labels.txt \
	--model-config example/models.json \
	--model        lr \
	--gene-sets    phase1-genesets.txt \
	--n-jobs       1 \
	--logdir       logs

# select candidate genes for each gene set
python bin/phase2-select.py \
	--gene-sets phase1-genesets.txt \
	--logdir    logs \
	--threshold 75 \
	--visualize \
	--outfile   phase2-genesets.txt

# select candidate genes using random forest
python bin/phase2-rf.py \
	--dataset   example_data.txt \
	--labels    example_labels.txt \
	--gene-sets phase1-genesets.txt \
	--n-jobs    1 \
	--threshold 75 \
	--visualize \
	--outfile   phase2-rf-genesets.txt
