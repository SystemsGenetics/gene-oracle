import numpy as np
import pandas as pd
import sys



def split_filename(filename):
	tokens = filename.split(".")

	return ".".join(tokens[:-1]), tokens[-1]



def load(filename):
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
		sys.exit(-1)



def save(filename, df):
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
		sys.exit(-1)
