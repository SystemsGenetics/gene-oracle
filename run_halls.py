#/usr/bin/python

# runs all the hallmark sets iteratively

import numpy as np
import os

#top_dir = os.system("cd /home/csheare/DeepGTEx/datasets/hallmark_subsets")

hall_subs = os.listdir("./datasets/hallmark_subsets")

for i in range(len(hall_subs)):
    os.system('./create-sets.py -d gtex -p ./datasets/hallmark_subsets/' + hall_subs[i] + ' -t 70 -r 30 ')
    #find features
    temp = np.fromfile('./datasets/hallmark_subsets/' + hall_subs[i] + '/Artery-Tibial/000_GTEX-111FC-0426-SM-5N9CV.dat', dtype=np.float32)
    n_features = temp.shape[0]
    #future:change arguments and get from command line
    os.system('python nn_gtex.py --n_input ' + str(n_features) + ' --n_classes 53 --batch_size 256 --lr 0.001 --epochs 1000')
