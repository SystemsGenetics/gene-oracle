#/usr/bin/python

# runs all the hallmark sets iteratively

import numpy as np
import os
import subprocess

#top_dir = os.system("cd /home/csheare/DeepGTEx/datasets/hallmark_subsets")

hall_subs = os.listdir("./datasets/hallmark_subsets")
hall_subs.sort()
files = ['random_accuracy.txt']
h1 = [1024]
h2 = [1024]
h3 = [1024]

accs = np.zeros((50,100))

for j in range(100):
	print('Iteration ' + str(j))
	out = []
	for i in range(len(hall_subs)):

		temp = np.fromfile('./datasets/hallmark_subsets/' + hall_subs[i] + '/Artery-Tibial/000_GTEX-111FC-0426-SM-5N9CV.dat', dtype=np.float32)
		n_features = temp.shape[0]

		os.system('python ./random_gene_generator.py --num ' + str(n_features))

		os.system('python ./create-sets.py -d gtex -p ./datasets/GTEx_Data_Random ' + ' -t 70 -r 30 ')
		#find features

		#future:change arguments and get from command line
		accs[i,j] = subprocess.check_output('python nn_gtex.py --n_input ' + str(n_features) + \
			' --n_classes 53 --batch_size 256 --lr 0.001 --epochs 75 --h1 ' + str(h1[0]) + ' --h2 ' + str(h2[0]) + ' --h3 ' + str(h3[0]), shell=True)


np.save('./logs/accuracies.npy', accs)

means = np.mean(accs, axis=1)

fp = open('./logs/' + files[0], 'w')
for val,name in zip(means,hall_subs):
    fp.write('{0}\t{1}'.format(name, val))    
fp.close() 
