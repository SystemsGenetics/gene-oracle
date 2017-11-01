#/usr/bin/python

import numpy as np 
import os

hallmarks = []
gtex_string_data = []
ensembl_map = []

print('loading hallmark gene list...')
with open('../BROAD_LISTS/h.all.v6.0.symbols.gmt') as f:
	content = f.readlines()

for c in content:
	hallmarks.append(np.array(c.split('\t')))

print('loading gtex GEM string data...')
gtex = np.load('../datasets/gtex_gct_data_string.npy')


subs = []
print('mapping genes...')

if not os.path.exists('../datasets/hallmark_numpys/'):
	os.makedirs('../datasets/hallmark_numpys/')

for h in hallmarks:
	subset = np.zeros((h.shape[0] - 2, gtex.shape[1]))
	subset = subset.astype(np.str)

	for c in range(2, h.shape[0]):
	
		a = np.where(gtex[:,1]==h[c])

		if (len(a[0]) != 0):
			subset[c - 2, :] = np.copy(gtex[a[0][0],:])
			
	np.save('../datasets/hallmark_numpys/' + h[0] + '.npy', subset)

