#/usr/bin/python

import numpy as np 
import os
import subprocess

sub = np.load('../datasets/hallmark_numpys/HALLMARK_HEDGEHOG_SIGNALING.npy')

genes = sub[:,1].tolist()

#top_dir = os.system("cd /home/csheare/DeepGTEx/datasets/hallmark_subsets")

files = ['hedgehog_single_gene_accuracy.txt']
h1 = [1024]
h2 = [1024]
h3 = [1024]

accs = np.zeros((35,1))

for j in range(30):
	for i in range(len(hall_subs)):
		os.system('python ../data_scripts/random_gene_generator.py --num ' + str(1))

		os.system('python ../data_scripts/create-sets.py -d gtex -p ../datasets/GTEx_Data_Random ' + ' -t 70 -r 30 ')
		#find features

		#future:change arguments and get from command line
		accs[i,j] = subprocess.check_output('python ../models/nn_gtex.py --n_input ' + str(n_features) + \
			' --n_classes 53 --batch_size 256 --lr 0.001 --epochs 75 --h1 ' + str(h1[0]) + ' --h2 ' + str(h2[0]) + ' --h3 ' + str(h3[0]), shell=True)


np.save('../logs/accuracies.npy', accs)

means = np.mean(accs, axis=1)

fp = open('../logs/' + files[0], 'w')
for val,name in zip(means,hall_subs):
    fp.write('{0}\t{1}'.format(name, val))    
fp.close() 
