# Gene Oracle

This repository contains the code for the Gene Oracle project. Gene Oracle is an ongoing research effort to discover biomarker genes using gene expression data. Gene Oracle identifies gene sets which provide the most predictive power, based on how well they classify samples in a gene expression dataset.

## Installation

All of Gene Oracle's dependencies can be installed via Anaconda. On a shared system (such as a university research cluster), it is recommended that you install everything in an Anaconda environment:

```bash
# specific to Clemson's Palmetto cluster
module add anaconda3/5.1.0

conda create -n gene-oracle python=3.6 tensorflow-gpu=1.13.1 matplotlib numpy pandas scikit-learn seaborn
```

You must then "activate" your environment in order to use it:
```bash
conda activate gene-oracle

# use gene-oracle

conda deactivate
```

After that, simply clone this repository to use Gene Oracle.
```bash
git clone https://github.com/SystemsGenetics/gene-oracle.git
cd gene-oracle

# run the example
example/run-example.sh
```

## Usage

Gene Oracle consists of two phases, (1) gene set analysis and (2) gene subset analysis. This process encompasses multiple scripts which are run in sequence. The easiest way to learn how to run these scripts, as well as the input / output data involved, is to run the example script as shown above. It demonstrates how to run Gene Oracle on synthetic input data from `make-input-data.py`.

### Input Data

Gene Oracle takes three primary inputs: (1) a gene expression matrix (GEM), (2) a list of sample labels, and (3) a list of gene sets. These inputs are described below.

The __gene expression matrix__ should be a plaintext file with rows being samples and columns being genes (features). Values in each row should be separated by tabs.
```
	Gene1	Gene2	Gene3	Gene4
Sample1	0.523	0.991	0.421	0.829
Sample2	8.891	7.673	3.333	9.103
Sample3	4.444	5.551	6.102	0.013
```

For large GEM files, it is recommended that you convert the GEM to numpy format using `convert.py` from the [GEMprep](https://github.com/SystemsGenetics/GEMprep) repo, as TSPG can load this binary format much more quickly than it does the plaintext format. The `convert.py` script can also transpose your GEM if it is arranged the wrong way:
```bash
python bin/convert.py GEM.txt GEM.npy --transpose
```

This example will create three files: `GEM.npy`, `GEM.rownames.txt`, and `GEM.colnames.txt`. The latter two files contain the row names and column names, respectively. Make sure that the rows are samples and the columns are genes!

The __label file__ should contain a label for each sample, corresponding to something such as a condition or phenotype state for the sample. This file should contain two columns, the first being the sample names and the second being the labels. Values in each row should be separated by tabs.
```
Sample1	Label1
Sample2	Label2
Sample3	Label3
Sample4	Label4
```

The __gene set list__ should contain the name and genes for a gene set on each line, similar to the GMT format. The gene names should be identical to those used in the GEM file. Values on each row should be separated by tabs.
```
GeneSet1	Gene1	Gene2	Gene3
GeneSet2	Gene2	Gene4	Gene5	Gene6
```

### Phase 1: Gene Set Analysis

The script `phase1-evaluate.py` takes a list of gene sets and evaluates each gene set by training and evaluating a classifier on the input dataset with only the genes in the set. This script can also evaluate the entire set of genes in the input dataset, as well as random gene sets.

The script `phase1-select.py` takes evaluation results for gene sets and compares them to results for random sets of equal size. It uses Welch's _t_-test (Student's _t_-test) to determine the statistical significance of a gene set's score as compared to a null distribution for the given set size. Larger gene sets tend to yield higher classification accuracies, so the _t_-test is used to eliminate this bias when selecting gene sets for subset analysis.

### Phase 2: Gene Subset Analysis

The script `phase2-evaluate.py` takes a list of gene sets and evaluates subsets of each gene set in order to determine the most salient genes in the gene set. This script can also analyze random gene sets in the same manner.

The script `phase2-select.py` takes evaluation results for the subsets selected by the previous script, measures the saliency of each gene by how frequently it appeared in all subsets, and separates "candidate" genes from "non-candidate" genes according to a threshold.

## Nextflow

### Dependencies

This repository also provides a Nextflow pipeline for running Gene Oracle. All you need is [nextflow](https://nextflow.io/), [Docker](https://docker.com/), and [nvidia-docker](https://github.com/NVIDIA/nvidia-docker). On HPC systems, you can use [Singularity](https://www.sylabs.io/singularity/) in lieu of Docker. If for some reason you can't use either container software, you will have to install Gene Oracle and its dependencies on your local machine.

### Input Data

The nextflow pipeline assumes you have your input data arranged as follows:
```
input/
  {dataset1}.emx.txt
  {dataset1}.labels.txt
  {dataset2}.emx.npy
  {dataset2}.rownames.txt
  {dataset2}.colnames.txt
  {dataset2}.labels.txt
  ...
  {genesets1}.genesets.txt
  {genesets2}.genesets.txt
  ...
```

This way, you can place as many gene subsets and datasets and the pipeline will process all of them in a single run.

### Usage

Here is a basic usage:
```
nextflow run systemsgenetics/gene-oracle
```

This example will download this pipeline to your machine and use the default `nextflow.config` in this repo. It will assume that you have Gene Oracle installed natively, and it will process all input files in the `input` directory, saving all output files to the `output` directory, as defined in `nextflow.config`.

You can also create your own `nextflow.config` file; nextflow will check for a config file in your current directory before defaulting to config file in this repo. You will most likely need to customize this config file as it provides options such as which experiments to run, how many chunks to use where applicable, and various other command-line parameters for Gene Oracle. The config file also allows you to define your own "profiles" for running this pipeline in different environments. Consult the Nextflow documentation for more information on what environments are supported.

To use Docker or Singularity, run nextflow with the `-with-docker` or `-with-singularity` flag. You can resume a failed run with the `-resume` flag. Consult the Nextflow documentation for more information on these and other options.

### Kubernetes

You can run this pipeline, as well as any other Nextflow pipeline, on a [Kubernetes](https://kubernetes.io/) cluster with minimal effort. Consult the [kube-runner](https://github.com/SystemsGenetics/kube-runner) repo for instructions.
