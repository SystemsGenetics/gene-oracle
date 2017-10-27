#/usr/bin/python

# runs all the hallmark sets iteratively

import numpy as np
import os
import subprocess

#top_dir = os.system("cd /home/csheare/DeepGTEx/datasets/hallmark_subsets")

hall_subs = os.listdir("../datasets/hallmark_subsets")
hall_subs.sort()
files = ['Xav_1024x1024x1024_100e.txt']
h1 = [1024]
h2 = [1024]
h3 = [1024]

for j in range(len(files)):
	out = []
	for i in range(len(hall_subs)):
	    os.system('../data_scripts/create-sets.py -d gtex -p ../datasets/hallmark_subsets/' + hall_subs[i] + ' -t 70 -r 30 ')
	    #find features
	    temp = np.fromfile('../data_scripts/datasets/hallmark_subsets/' + hall_subs[i] + '/Artery-Tibial/000_GTEX-111FC-0426-SM-5N9CV.dat', dtype=np.float32)
	    n_features = temp.shape[0]
	    #future:change arguments and get from command line
	    out.append(subprocess.check_output('python ../models/nn_gtex_reinforce.py --n_input ' + str(n_features) + \
	    	' --n_classes 53 --batch_size 16 --lr 0.001 --epochs 100 --h1 ' + str(h1[j]) + ' --h2 ' + str(h2[j]) + ' --h3 ' + str(h3[j]), shell=True))

	fp = open('../logs/' + files[j], 'w')
	for val,name in zip(out,hall_subs):
	    fp.write('{0}\t{1}'.format(name, val))    
	fp.close() 
