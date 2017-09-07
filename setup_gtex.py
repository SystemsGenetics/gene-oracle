#/usr/bin/python

# class to handle gtex data for NN training/testing
import os
import numpy as np
from random import randint
from shutil import copyfile

def get_subdirs(dir):
    return [name for name in os.listdir(dir)
            if os.path.isdir(os.path.join(dir, name))]

def get_files(path):
	return [name for name in os.listdir(path)
		if os.path.isfile(os.path.join(path, name))]

def balance_dataset(num_samples, dataset_path):
	subdirs = get_subdirs(dataset_path)
	subdirs.sort()
	classes = subdirs

	for c in classes:
		path = dataset_path + '/' + c
		files = get_files(path)

		if len(files) > num_samples:
			while len(get_files(path)) > num_samples:
				os.remove(path + '/' + files.pop())
		elif len(files) < num_samples:
			i = 0
			while len(get_files(path)) < num_samples:
				files = get_files(path)
				f = files[randint(0, len(files) - 1)]
				copyfile(path + '/' + f, path + '/' + str(i) + f)
				i = i + 1
		else:
			continue

class train_t(object):
	def __init__(self, labels, data):
		self.labels = labels
		self.data = np.transpose(data)
		self.num_examples = data.shape[1]

	def next_batch(self, batch_size, index):
		idx = index * batch_size
		n_idx = index * batch_size + batch_size
		return self.data[idx:n_idx, :], self.labels[idx:n_idx, :]

class test_t(object):
	def __init__(self, labels, data):
		self.labels = labels
		self.data = np.transpose(data)
		self.num_examples = data.shape[1]

	def next_batch(self, batch_size, index):
		idx = index * batch_size
		n_idx = index * batch_size + batch_size
		return self.data[idx:n_idx, :], self.labels[idx:n_idx, :]

class GTEx(object):
	def __init__(self, dataset_path, train_path, test_path):
		self.dataset_path = dataset_path
		self.train_path = train_path
		self.test_path = test_path
		self.classes = self.get_classes(dataset_path)
		self.train = self.get_train(train_path)
		self.test = self.get_test(test_path)

	def get_train(self, train_path):
		files = get_files(train_path)	
		labels = np.zeros((len(files), len(self.classes)))

		i = 0
		data = np.zeros((56238, len(files)))

		for file in files:
			tmp = np.fromfile(os.path.join(train_path, file), dtype=np.float32)
			data[:,i] = np.copy(tmp)
			i = i + 1

		i = 0
		for i in range(len(files)):
			tmp = files[i].split('_')
			labels[i, self.classes.index(tmp[0])] = 1

		train = train_t(labels, data)

		return train

	def get_test(self, test_path):
		files = get_files(test_path)	
		labels = np.zeros((len(files), len(self.classes)))

		data = np.zeros((56238, len(files)))
		i = 0
		for file in files:
			tmp = np.fromfile(os.path.join(test_path, file), dtype=np.float32)
			data[:,i] = np.copy(tmp)
			i = i + 1

		i = 0
		for i in range(len(files)):
			tmp = files[i].split('_')
			labels[i, self.classes.index(tmp[0])] = 1

		test = test_t(labels, data)

		return test

	def get_classes(self, dataset_path):
		subs = get_subdirs(dataset_path)
		subs.sort()
		return subs
