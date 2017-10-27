#/usr/bin/python


import numpy as np
from sklearn import preprocessing
import matplotlib.pyplot as plt
import os

if os.path.exists('../datasets/gtex.npy'):
	data = np.load('../datasets/gtex.npy')
else:
	classes = os.listdir('../datasets/hallmark_subsets/HALLMARK_ADIPOGENESIS')

	data = np.zeros((200,8555))
	i = 0

	for c in classes:
		files = os.listdir('../datasets/hallmark_subsets/HALLMARK_ADIPOGENESIS/' + c)
		for f in files:
		    sample = np.fromfile('../datasets/hallmark_subsets/HALLMARK_ADIPOGENESIS/' + c + '/' + f, dtype=np.float32)
		    data[:, i] = sample
		    i = i + 1

	np.save('../datasets/gtex.npy', data)

#scale data
maxabs = preprocessing.MaxAbsScaler()
out = maxabs.fit_transform(data)

# plot graph
fig, ax = plt.subplots()
heatmap = ax.pcolor(out, cmap='coolwarm')
#plt.colorbar()
plt.show()
