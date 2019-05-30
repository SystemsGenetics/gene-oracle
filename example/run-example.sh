#!/bin/bash
# Example usage of gene-oracle on a synthetic dataset.

# create synthetic input data
python scripts/make-classification.py

# evaluate gene sets
python scripts/phase1-evaluate.py \
	--dataset      example_data.txt \
	--labels       example_labels.txt \
	--model-config example/models.json \
	--gene-sets    example_genesets.txt \
	--num-folds    5 \
	--outfile      phase1-curated.txt

python scripts/phase1-evaluate.py \
	--dataset      example_data.txt \
	--labels       example_labels.txt \
	--model-config example/models.json \
	--random \
	--random-range 1 20 1 \
	--random-iters 10 \
	--num-folds    5 \
	--outfile      phase1-random.txt

# perform Welch's t-test on gene sets
python scripts/phase1-select.py \
	--scores-fg phase1-curated.txt \
	--scores-bg phase1-random.txt \
	--gene-sets example_genesets.txt

# perform combinatorial analysis on gene sets
python scripts/phase2-evaluate.py \
	--dataset      example_data.txt \
	--labels       example_labels.txt \
	--model-config example/models.json \
	--gene-sets    example_genesets.txt \
	--random \
	--logdir       logs

# select candidate genes for each gene set
python scripts/phase2-select.py \
	--gene-sets example_genesets.txt \
	--logdir    logs
