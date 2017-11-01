# DeepGTEx
This repository contains all code, scripts, and other for the DeepGTEx project. DeepGTEx is an ongoing research effort to effectively classify RNA expression level data from the Genotype-Tissue Expression project, located [here](https://gtexportal.org/home/). Running models from this repo requires a downloaded version of the RNA-seq data Gene TPM's, located [here](https://gtexportal.org/home/datasets) (a sign-in is required to download the dataset for privacy reasons). 






    usage: nn_gtex.py [-h] [--lr LR] [--epochs EPOCHS] [--h1 H1] [--h2 H2]
                      [--h3 H3] [--batch_size BATCH_SIZE]
                      [--display_step DISPLAY_STEP] [--n_input N_INPUT]
                      [--n_classes N_CLASSES] [--beta BETA] [--load LOAD]
                      [--confusion CONFUSION]

    optional arguments:
      -h, --help              show this help message and exit
      --lr FLOAT              learning rate
      --epochs INT            no. of training epoch
      --h1 H1                 no. of neurons in hidden layer 1
      --h2 H2                 no. of neurons in hidden layer 2
      --h3 H3                 no. of neurons in hidden layer 3
      --batch_size INT        batch size
      --display_step INT      print updates after this many steps
      --n_input INT           number of input features
      --n_classes INT         number of classes
      --beta BETA             hyperparemeter for l1 regularization of weights
      --load LOAD             load weights from previous run
      --confusion BOOL        generate confusion matrix (1) or no (0)
  
