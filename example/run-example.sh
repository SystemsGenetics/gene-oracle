#!/bin/bash
# Example usage of gene-oracle on a synthetic dataset.

# create synthetic input data
python scripts/make-classification.py

# evaluate gene sets
python scripts/phase1-evaluate.py \
	--dataset      example_data.txt \
	--labels       example_labels.txt \
	--model_config example/models.json \
	--gene_sets    example_genesets.txt \
	--num_folds    5 \
	--outfile      phase1-genesets.txt

python scripts/phase1-evaluate.py \
	--dataset      example_data.txt \
	--labels       example_labels.txt \
	--model_config example/models.json \
	--random \
	--random_range 1 20 1 \
	--random_iters 10 \
	--num_folds    5 \
	--outfile      phase1-random.txt

# perform Welch's t-test on gene sets
python scripts/phase1-select.py \
	--random    phase1-random.txt \
	--subset    phase1-genesets.txt \
	--gene_sets example_genesets.txt

# perform combinatorial analysis on gene sets
python scripts/phase2-evaluate.py \
	--dataset      example_data.txt \
	--labels       example_labels.txt \
	--model_config example/models.json \
	--gene_sets    example_genesets.txt \
	--random \
	--logdir       logs
