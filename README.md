# DeepGTEx
This repository contains all code, scripts, and other for the DeepGTEx project. DeepGTEx is an ongoing research effort to effectively classify RNA expression level data from the Genotype-Tissue Expression project, located [here](https://gtexportal.org/home/). Running models from this repo requires a downloaded version of the RNA-seq data Gene TPM's, located [here](https://gtexportal.org/home/datasets) (a sign-in is required to download the dataset for privacy reasons). 

## Requirements and Setup
Create a virtual environment with anaconda3 on Clemson's Palmetto Cluster:
    
    module add anaconda/4.3.0
    conda create -n {name of env} python=2.7
    
To activate your environment, simply do:

    source activate {name of env}
    
Once a virtualenv is created, the following software is required to run the models:

    conda install -n yourenvname [package]
    example: conda install -n gtex tensorflow-gpu==1.3.0
    
    tensorflow-gpu (1.3.0)
    scikit-learn (0.19.0)
    numpy (1.13.1)
    argparse (1.4.0)
    
## Usage
classify.py can be used to classify a dataset using all the features, a random subset of features, or a specified subset of features.  
Users are required to input a path to three files:
*a numpy data array that has features (genes) row-wise and samples column wise
*a numpy data array that contains a list of every gene (str) in the exact order as the dataset
*a json file that contains the number of samples per class
    *note: the dataset is assumed to be in order of the json file

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

  
## Feature Engineering
To determine the relevance of a particular gene or a subgroup of genes, [subset_gene_test.py](https://github.com/CUFCTL/DeepGTEx/blob/master/scripts/subset_gene_test.py) can be used to generate subgroups from hallmark sets. 
