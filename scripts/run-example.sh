#!/bin/bash
# Example usage of gene-oracle on a synthetic dataset.

set -e

DATASET="example.emx.txt"
LABELS="example.labels.txt"
GMT_FILE="example.genesets.txt"
MODEL_CONFIG="models.json"
MODEL="lr"
OUTPUT_DIR="output"

# initialize environment
source activate gene-oracle

# initialize output directory
rm -rf ${OUTPUT_DIR}

mkdir -p ${OUTPUT_DIR}

# create synthetic input data
python bin/make-inputs.py \
    --n-samples 1000 \
    --n-genes   100 \
    --n-classes 5 \
    --n-sets    10 \
    --visualize

# evaluate gene sets
python bin/phase1-evaluate.py \
    --dataset      ${DATASET} \
    --labels       ${LABELS} \
    --gene-sets    ${GMT_FILE} \
    --model-config ${MODEL_CONFIG} \
    --model        ${MODEL} \
    --random \
    --random-iters 10 \
    --cv           5 \
    --n-jobs       1 \
    --visualize \
    --output-dir   ${OUTPUT_DIR}

# select gene sets which score higher over random
python bin/phase1-select.py \
    --dataset    ${DATASET} \
    --gene-sets  ${GMT_FILE} \
    --scores     ${OUTPUT_DIR}/phase1-scores.txt \
    --threshold  1 \
    --n-sets     5 \
    --visualize \
    --output-dir ${OUTPUT_DIR}

# perform combinatorial analysis on selected gene sets
# python bin/phase2-evaluate.py \
#     --dataset      ${DATASET} \
#     --labels       ${LABELS} \
#     --gene-sets    ${OUTPUT_DIR}/phase1-genesets.txt \
#     --model-config ${MODEL_CONFIG} \
#     --model        ${MODEL} \
#     --n-jobs       1 \
#     --logdir       ${OUTPUT_DIR}/logs

# select candidate genes for each gene set
# python bin/phase2-select.py \
#     --gene-sets  ${OUTPUT_DIR}/phase1-genesets.txt \
#     --logdir     ${OUTPUT_DIR}/logs \
#     --threshold  75 \
#     --visualize  \
#     --output-dir ${OUTPUT_DIR}

# select candidate genes using random forest
python bin/phase2-rf.py \
    --dataset    ${DATASET} \
    --labels     ${LABELS} \
    --gene-sets  ${OUTPUT_DIR}/phase1-genesets.txt \
    --n-jobs     1 \
    --threshold  75 \
    --visualize  \
    --output-dir ${OUTPUT_DIR}
