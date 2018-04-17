#/usr/bin/python

# autoencoder implementation for compressing sets of gene data

import numpy as np 
import tensorflow as tf 
from sklearn import preprocessing
import sys, os, argparse

sys.path.append(os.path.dirname(os.getcwd()))


def next_batch(data, batch_size, index):
	idx = index * batch_size
	n_idx = index * batch_size + batch_size
	return data[idx:n_idx, :]


class autoencoder:

	def __init__ (self, epochs=100, lr=0.01, n_input=56202, batch_size=16):
		self.epochs = epochs
		self.lr = lr
		self.n_input = n_input
		self.batch_size = batch_size


	def encoder(self, x, weights, biases):
		l1 = tf.nn.relu(tf.add(tf.matmul(x, weights['encoder_h1']), \
											biases['encoder_b1']))
		l2 = tf.nn.relu(tf.add(tf.matmul(l1, weights['encoder_h2']), \
											biases['encoder_b2']))
		return l2


	def decoder(self, x, weights, biases):
		l1 = tf.nn.relu(tf.add(tf.matmul(x, weights['decoder_h1']), \
											biases['decoder_b1']))
		l2 = tf.nn.relu(tf.add(tf.matmul(l1, weights['decoder_h2']), \
											biases['decoder_b2']))
		return l2



	def run(self, data):
		tf.reset_default_graph()

		x = tf.placeholder("float", [None, self.n_input])

		weights = {
			'encoder_h1': tf.get_variable('encoder_h1', shape=[self.n_input, 1000],  initializer=tf.keras.initializers.glorot_uniform()),
			'encoder_h2': tf.get_variable('encoder_h2', shape=[1000, 200], initializer=tf.keras.initializers.glorot_uniform()),
			'decoder_h1': tf.get_variable('decoder_h1', shape=[200, 1000], initializer=tf.keras.initializers.glorot_uniform()),
			'decoder_h2': tf.get_variable('decoder_h2', shape=[1000, self.n_input], initializer=tf.keras.initializers.glorot_uniform())
		}

		biases = {
			'encoder_b1': tf.get_variable('encoder_b1', shape=[1000], initializer=tf.keras.initializers.glorot_uniform()),
			'encoder_b2': tf.get_variable('encoder_b2', shape=[200], initializer=tf.keras.initializers.glorot_uniform()),
			'decoder_b1': tf.get_variable('decoder_b1', shape=[1000], initializer=tf.keras.initializers.glorot_uniform()),
			'decoder_b2': tf.get_variable('decoder_b2', shape=[self.n_input], initializer=tf.keras.initializers.glorot_uniform())
		}

		# run architecture
		encode = self.encoder(x, weights, biases)
		decode = self.decoder(encode, weights, biases)

		# set up output/label
		y_pred = decode
		y_true = x

		# define loss function and optimizer
		loss = tf.reduce_mean(tf.pow(y_true - y_pred, 2))
		optimizer = tf.train.RMSPropOptimizer(self.lr).minimize(loss)

		# Initialize the variables
		init = tf.global_variables_initializer()

		# Add ops to save and restore all the variables.
		saver = tf.train.Saver()



		# create a summary for our cost and accuracy
		tf.summary.scalar("cost", loss)

		# start tensorflow sess
		with tf.Session() as sess:

			sess.run(init)

			# set up writer for tensorboard
			writer = tf.summary.FileWriter('./tf_log', sess.graph)

			for epoch in range(1, self.epochs + 1):

				avg_cost = 0.0
				total_batch = int(data.shape[0]/self.batch_size)

				for i in range(total_batch):

					batch_x = next_batch(data, self.batch_size, i)

					_, l = sess.run([optimizer, loss], feed_dict={x: batch_x})

					#print ('batch ' + str(i) + '/' + str(total_batch) + ' loss: ' + str(l))

					# Compute average loss
					avg_cost += l / total_batch

				print("Epoch: ", '%04d' % (epoch), "loss: ", "{:.9f}".format(avg_cost))

			saver.save(sess, 'weights/test_autoencoder.ckpt')

























