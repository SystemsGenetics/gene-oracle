#/usr/bin/python

# runs all the hallmark sets iteratively

import numpy as np
import os
import subprocess

#top_dir = os.system("cd /home/csheare/DeepGTEx/datasets/hallmark_subsets")
batch_size = [512, 512, 512, 512, 256, 256, 256, 128, 128, 64, 64, 32, 32]
datasets   = [5, 15, 25, 50, 75, 100, 150, 200, 250, 350, 500, 750, 1000] #1500, 2500, 5000, 10000, 56238]
files = ['r_accuracy_varying_features.txt']
h1 = [1024]
h2 = [1024]
h3 = [1024]

accs = np.zeros((13,30))

for j in range(30):
	print('Iteration ' + str(j))
	out = []
	for i in range(len(datasets)):
		print('\tIter ' + str(i))

		# temp = np.fromfile('./datasets/hallmark_subsets/' + hall_subs[i] + '/Artery-Tibial/000_GTEX-111FC-0426-SM-5N9CV.dat', dtype=np.float32)
		# n_features = temp.shape[0]

		os.system('python ./random_gene_generator.py --num ' + str(datasets[i]))

		os.system('python ./create-sets.py -d gtex -p ../datasets/GTEx_Data_Random ' + ' -t 70 -r 30 ')
		#find features

		accs[i,j] = subprocess.check_output('python ../models/nn_gtex.py --n_input ' + str(datasets[i]) + \
				' --n_classes 53 --batch_size ' + str(batch_size[i]) + ' --lr 0.001 --epochs 75 --h1 ' + str(h1[0]) + ' --h2 ' + str(h2[0]) + ' --h3 ' + str(h3[0]), shell=True)

	np.save('../logs/accuracies_rand.npy', accs)

means = np.mean(accs, axis=1)

fp = open('../logs/' + files[0], 'w')
for val,name in zip(means,datasets):
    fp.write('{0}\t{1}\n'.format(name, val))    
fp.close() 
