#!/bin/bash
# Example usage of gene-oracle on a synthetic dataset.

# create synthetic input data
python scripts/make-classification.py

# evaluate gene sets
python scripts/phase1-evaluate.py \
	--dataset      example_data.txt \
	--labels       example_labels.txt \
	--model-config example/models.json \
	--model        mlp-tf \
	--gene-sets    example_genesets.txt \
	--random \
	--random-iters 10 \
	--cv           5 \
	--n-jobs       1 \
	--outfile      phase1-scores.txt

# select gene sets which score higher over random
python scripts/phase1-select.py \
	--scores    phase1-scores.txt \
	--gene-sets example_genesets.txt

# perform combinatorial analysis on selected gene sets
python scripts/phase2-evaluate.py \
	--dataset      example_data.txt \
	--labels       example_labels.txt \
	--model-config example/models.json \
	--model        mlp-tf \
	--gene-sets    example_genesets.txt \
	--logdir       logs

# select candidate genes for each gene set
python scripts/phase2-select.py \
	--gene-sets example_genesets.txt \
	--logdir    logs
