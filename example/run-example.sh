#!/bin/bash

# create synthetic input data
python scripts/make-classification.py

# evaluate gene sets
python scripts/phase1-evaluate.py \
	--dataset      example_data.txt \
	--labels       example_labels.txt \
	--model_config example/model_config.json \
	--gene_sets    example_genesets.txt \
	--num_folds    5 \
	--outfile      phase1-genesets.txt

python scripts/phase1-evaluate.py \
	--dataset      example_data.txt \
	--labels       example_labels.txt \
	--model_config example/model_config.json \
	--random \
	--random_range 1 20 \
	--random_iters 10 \
	--num_folds    5 \
	--outfile      phase1-random.txt

# perform Welch's t-test on gene sets
python scripts/phase1-screen.py \
	--random    phase1-random.txt \
	--subset    phase1-genesets.txt \
	--gene_sets example_genesets.txt

# perform candidate selection on gene sets
# python scripts/gene-oracle.py # ...
