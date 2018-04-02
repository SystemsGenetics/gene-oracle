#
# FILE: GTEx.py
# USE:  Creates object that holds data and labels for gene specified dataset
#

import numpy as np
import random
import sys

class data_t(object):
	def __init__(self, data, labels):
		self.labels = labels
		self.data = data
		self.num_examples = data.shape[0]

	def next_batch(self, batch_size, index):
		idx = index * batch_size
		n_idx = index * batch_size + batch_size
		return self.data[idx:n_idx, :], self.labels[idx:n_idx, :]

class GTEx:
	def __init__(self, data, total_gene_list=None, sub_gene_list=None, train_split=70, test_split=30):
		self.num_classes = len(data)
		self.train, self.test = self.split_set(data, total_gene_list, sub_gene_list, train_split, test_split)


	# 
	# USAGE:
	# 		create a new data dictionary with classes as keys that contain a subsample of original data,
	#		specified by the 'sub_gene_list' 
	# PARAMS:
	#		orig_data:       dictionary containing classes as keys, with values as matrix of class samples
	#						 with all possible genes 
	#		total_gene_list: list of every gene in the original data GEM
	#   	sub_gene_list:   specified genes from a subset of total gene list
	#
	def extract_requested_genes(self, orig_data, total_gene_list, sub_gene_list):
		
		#get genes in sub_gen_list from total_gene_list
		gene_indexes = []
		for i in range(len(sub_gene_list)):
			gene_indexes.append(np.argwhere(total_gene_list == sub_gene_list[i]))

		# dictionary for requested data
		req_data = {}

		# iterate through dictionary, replace old data matrix with reduced data matrix
		for k in sorted(orig_data.keys()):
			reduced_data = np.zeros((len(gene_indexes), orig_data[k].shape[1]))

			for idx in xrange(0, len(gene_indexes)):
				reduced_data[idx] = orig_data[k][gene_indexes[idx][0],:]

			req_data[k] = reduced_data

		return req_data


	#
	# USAGE:
	#	shuffle the data and labels in the same order for a data set and transform the data and labels into numpy arrays
	# PARAMS:
	#	data:	the data values for a dataset
	# 	labels: the labels associated with data for a dataset
	#
	def shuffle_and_transform(self, data,labels):
		new_data = []
		new_labels = []
		
		samples = random.sample(xrange(len(data)),len(data))
		
		for i in samples:
			new_data.append(data[i])
			new_labels.append(labels[i])

		# convert lists to numpy arrays
		np_data = np.asarray(new_data)
		np_labels = np.asarray(new_labels)
		
		return data_t(np_data,np_labels)



	#
	# USAGE:
	#       split the data in data dictionary to the specified weights
	# PARAMS:
	#		data: 			 dictionary containing all samples of each class, with the class name being the key
	#       total_gene_list: list of every sorted gene in the GEM
	#       sub_gene_list:   list of genes in the subset wished to be extracted. this can be the total gene list
	#						 if all genes are wanting to be used
	#       train_split:     percentage (in integer form) for training class split
	#       test_split:      percentage (in integer form) for testing class split
	#
	def split_set(self, data, total_gene_list=None, sub_gene_list=None, train_split=70, test_split=30):

		if test_split + train_split != 100:
			print('Test and train split must sum to 100!')
			sys.exit(1)

		train_data = []
		train_labels = []
		test_data = []
		test_labels = []

		if sub_gene_list is not None:
			data = self.extract_requested_genes(data, total_gene_list, sub_gene_list)

		idx = 0 # keep count of index in dictionary

		# gather training and testing examples and labels into a list by randomly selecting indices 
		# of the amount of data in each class
		for k in sorted(data.keys()):
			num_train = int(data[k].shape[1] * train_split / 100)

			samples = random.sample(xrange(data[k].shape[1]),data[k].shape[1])
			samples_train = samples[0:num_train]
			samples_test = samples[num_train:]

			for i in xrange(len(samples_train)):
				train_data.append(data[k][:,samples_train[i]])

				label = np.zeros(self.num_classes)
				label[idx] = 1

				train_labels.append(label)

			for i in xrange(len(samples_test)):
				test_data.append(data[k][:,samples_test[i]])

				label = np.zeros(self.num_classes)
				label[idx] = 1

				test_labels.append(label)

			idx = idx + 1

		train = self.shuffle_and_transform(train_data, train_labels)
		test = self.shuffle_and_transform(test_data, test_labels)

		return [train, test]

