import copy
import json
import numpy as np
import pandas as pd
import sklearn.dummy
import sklearn.ensemble
import sklearn.linear_model
import sklearn.model_selection
import sklearn.neighbors
import sklearn.neural_network
import sklearn.preprocessing
import sklearn.svm
import sys

import models



def split_filename(filename):
	tokens = filename.split(".")

	return ".".join(tokens[:-1]), tokens[-1]



def load_dataframe(filename):
	basename, ext = split_filename(filename)

	if ext == "txt":
		# load dataframe from plaintext file
		return pd.read_csv(filename, index_col=0, sep="\t")
	elif ext == "npy":
		# load data matrix from binary file
		X = np.load(filename)

		# load row names and column names from text files
		rownames = np.loadtxt("%s_rownames.txt" % basename, dtype=str)
		colnames = np.loadtxt("%s_colnames.txt" % basename, dtype=str)

		# combine data, row names, and column names into dataframe
		return pd.DataFrame(X, index=rownames, columns=colnames)
	else:
		print("error: filename %s is invalid" % (filename))
		sys.exit(1)



def save_dataframe(filename, df):
	basename, ext = split_filename(filename)

	if ext == "txt":
		# save dataframe to plaintext file
		df.to_csv(filename, sep="\t", na_rep="NA", float_format="%.8f")
	elif ext == "npy":
		# save data matrix to binary file
		np.save(filename, np.array(df.values, dtype=np.float32, order="F"))

		# save row names and column names to text files
		np.savetxt("%s_rownames.txt" % basename, df.index, fmt="%s")
		np.savetxt("%s_colnames.txt" % basename, df.columns, fmt="%s")
	else:
		print("error: filename %s is invalid" % (filename))
		sys.exit(1)



def load_labels(filename):
	# load labels file
	labels = pd.read_csv(filename, sep="\t", header=None, index_col=0)

	# convert categorical labels to numerical labels
	labels = labels[1].values
	labels = sklearn.preprocessing.LabelEncoder().fit_transform(labels)

	return labels



def load_gene_sets(filename):
	# load file into list
	lines = [line.strip() for line in open(filename, "r")]
	lines = [line.split("\t") for line in lines]

	# map each gene set into a tuple of the name and genes in the set
	gene_sets = [(line[0], line[1:]) for line in lines]

	return gene_sets



def load_classifier(config_file, name):
	# load model params
	config = json.load(open(config_file))
	params = config[name] if name in config else {}

	# define dictionary of classifiers
	classifiers = {
		"dummy": sklearn.dummy.DummyClassifier,
		"knn": sklearn.neighbors.KNeighborsClassifier,
		"lr": sklearn.linear_model.LogisticRegression,
		"mlp": sklearn.neural_network.MLPClassifier,
		"mlp-tf": models.TensorflowMLP,
		"rf": sklearn.ensemble.RandomForestClassifier,
		"svm": sklearn.svm.SVC
	}

	# initialize classifier
	if name in classifiers:
		clf = classifiers[name](**params)

	# or print error if model name is invalid
	else:
		print("error: classifier \"%s\" not recognized" % name)
		sys.exit(1)

	return clf



def evaluate_gene_set(data, labels, clf, genes, cv=None, n_jobs=None):
	# extract dataset
	X = data[genes]

	# normalize dataset
	X = sklearn.preprocessing.MaxAbsScaler().fit_transform(X)

	# perform a single train/test split if cv is not specified
	if cv == None or cv == 1:
		# create train/test sets
		X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(X, labels, test_size=0.3)

		# evaluate gene set
		clf = copy.deepcopy(clf)
		clf.fit(X_train, y_train)

		return [clf.score(X_test, y_test)]

	# otherwise use cross-validation
	else:
		# evaluate gene set
		scores = sklearn.model_selection.cross_val_score(clf, X, y=labels, cv=cv, n_jobs=n_jobs)

		return list(scores)
