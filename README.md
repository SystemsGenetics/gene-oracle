# Gene Oracle

This repository contains all code, scripts, and other for the Gene Oracle project. Gene Oracle is an ongoing research effort to discover biomarker genes given RNA expression data. The data is expected to be preprocessed and in numpy float format with rows being genes (features) and columns being samples (each data point being an RNA expression level). Additionally, a list of genes, also in numpy format, is expected as input in the same order in which the float data is constructed. Finally, a json file containing class names and number of samples per class, in sample order, is required. 

## Requirements and Setup

Create a virtual environment with anaconda3 (this example assumes Clemson's Palmetto Cluster, but any compute cluster should work):
```bash
module add anaconda3/5.1.0
conda create -n gene-oracle python=3.5 tensorflow-gpu=1.8.0 matplotlib numpy pandas scikit-learn
```

To activate your environment, simply do:
```bash
conda activate gene-oracle
```

To deactivate your environment, simply do:
```bash
conda deactivate
```

## Usage

The script `classify.py` can be used to classify a dataset using all the features, a random subset of features, or a specified subset of features.  

Users are required to input a path to three files:
* a numpy data array that has features (genes) row-wise and samples column wise
* a numpy data array that contains a list of every gene (str) in the exact order as the dataset
* a json file that contains the number of samples per class
    * note: the dataset is assumed to be in order of the json file
    * Example: {"Adipose-Subcutaneous": 350, "Adipose-Visceral": 227,...,"Whole-Blood": 393}
    so for class Adipose-Subcutaneous there are 350 samples etc etc
* if a subset list is input, it is expected to be of the following format:
```
SetName1,Gene1,Gene2,Gene3
SetName2,Gene2,Gene4,Gene5,Gene6
...
```
* note: working towards a more robust data loading scheme

The following contains the example usage:
```
usage: classify.py [-h] --dataset DATASET --gene_list GENE_LIST --sample_json
                   SAMPLE_JSON --config CONFIG --out_file OUT_FILE
                   [--subset_list SUBSET_LIST] [--random_test]
                   [--num_random_genes NUM_RANDOM_GENES [NUM_RANDOM_GENES ...]]
                   [--rand_iters [RAND_ITERS]]

arguments:
  -h, --help            show this help message and exit
  --dataset DATASET     dataset to be used
  --gene_list GENE_LIST
                        list of genes in dataset (same order as dataset)
  --sample_json SAMPLE_JSON
                        json file containing number of samples per class
  --config CONFIG       json file containing network specifications
  --out_file OUT_FILE   output file to send results to
  --subset_list SUBSET_LIST
                        gmt/gct file containing subsets
  --random_test         Perform random test
  --num_random_genes NUM_RANDOM_GENES [NUM_RANDOM_GENES ...]
                        Number of random genes to assess
  --rand_iters [RAND_ITERS]
                        Number of iterations to perform for random
                        classification
```

## Candidate Selection

To determine the relevance of a particular gene or a subgroup of genes, `gene-oracle.py` can be used to generate subgroups from sets of genes. This file expects data to be the in the same format as the classification script above. 
```
usage: gene-oracle.py [-h] --dataset DATASET --gene_list GENE_LIST
                      --sample_json SAMPLE_JSON --config CONFIG
                      [--subset_list SUBSET_LIST] --set SET
                      [--num_genes NUM_GENES] --log_dir LOG_DIR

Run tests on specified subsets of a hallmark or random set

optional arguments:
  -h, --help            show this help message and exit
  --dataset DATASET     dataset to be used
  --gene_list GENE_LIST
                        list of genes in dataset (same order as dataset)
  --sample_json SAMPLE_JSON
                        json file containing number of samples per class
  --config CONFIG       json file containing network specifications
  --subset_list SUBSET_LIST
                        gmt/gct/txt file containing subsets
  --set SET             subset to be used
  --num_genes NUM_GENES
                        number of genes
  --log_dir LOG_DIR     directory where logs are stored
```
The file expects a set to be the name of group of genes in the subset_list file. This is the set that is decomposed to subsets according to phase 2 of the algorithm. Additionally, input the number of genes in the set and a log directory to put the output files in.
