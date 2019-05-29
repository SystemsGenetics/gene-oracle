import numpy as np
import os
import sklearn.base
import sys
import tensorflow as tf



class MLP(sklearn.base.BaseEstimator):

	def __init__(self, layers=[512,512,512], activations=["relu", "relu", "relu"], \
		dropout=False, lr=0.001, epochs=75, batch_size=16, \
		load=False, save=False, verbose=False):

		self.layers = layers
		self.activations = activations
		self.dropout = dropout
		self.lr = lr
		self.epochs = epochs
		self.batch_size = batch_size
		self.load = load
		self.save = save
		self.verbose = verbose



	def __del__(self):
		if hasattr(self, "_session"):
			self._session.close()


	def _initialize(self, n_input, n_classes):
		os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"

		tf.reset_default_graph()

		x = tf.placeholder("float", [None, n_input])
		y = tf.placeholder("float", [None, n_classes])

		# initialize weights and biases
		layer_sizes = [n_input] + [n for n in self.layers] + [n_classes]
		weights = {}
		biases = {}

		for i in range(1, len(self.layers) + 1):
			w = "h%d" % i
			b = "b%d" % i
			weights[w] = tf.get_variable(w, shape=[layer_sizes[i - 1], layer_sizes[i]], initializer=tf.contrib.layers.xavier_initializer())
			biases[b] = tf.get_variable(b, shape=[layer_sizes[i]], initializer=tf.contrib.layers.xavier_initializer())

		weights["out"] = tf.get_variable("out_w", shape=[self.layers[-1], n_classes], initializer=tf.contrib.layers.xavier_initializer())
		biases["out"] = tf.get_variable("out_b", shape=[n_classes], initializer=tf.contrib.layers.xavier_initializer())

		# initialize computational graph
		model = x

		for i in range(1, len(self.layers) + 1):
			# append weights and biases
			w = "h%d" % i
			b = "b%d" % i
			model = tf.add(tf.matmul(model, weights[w]), biases[b])

			# append activation function
			activation = self.activations[i - 1]

			if activation == "relu":
				model = tf.nn.relu(model)
			elif activation == "sigmoid":
				model = tf.nn.sigmoid(model)
			else:
				print("error: unrecognized activation function %s" % activation)
				sys.exit(1)

			# append dropout layer
			if self.dropout:
				model = tf.nn.dropout(model, 0.5)

		# append output layer
		logits = tf.add(tf.matmul(model, weights["out"]), biases["out"])
		model = tf.nn.softmax(logits)

		# define learning rate, loss, optimizer
		global_step = tf.Variable(0, trainable=False)
		learning_rate = tf.train.exponential_decay(self.lr, global_step=global_step, decay_steps=500, decay_rate=0.96, staircase=True)
		loss = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits_v2(logits=logits, labels=y))
		optimizer = tf.train.AdamOptimizer(learning_rate=learning_rate).minimize(loss, global_step=global_step)

		# save tf objects
		self._model = model
		self._lr = learning_rate
		self._loss = loss
		self._optimizer = optimizer
		self._x = x
		self._y = y

		# initialize computational graph
		if not hasattr(self, "_session"):
			self._session = tf.Session()

		init = tf.global_variables_initializer()
		self._session.run(init)



	def _onehot_encode(self, y):
		return np.array([self._classes == y_i for y_i in y])



	def _shuffle(self, x, y):
		indices = np.arange(x.shape[0])
		np.random.shuffle(indices)

		return x[indices], y[indices]



	def _next_batch(self, x, y, batch_size, index):
		a = index * batch_size
		b = index * batch_size + batch_size

		return x[a:b], y[a:b]



	def fit(self, x, y):
		n_samples = x.shape[0]
		n_input = x.shape[1]

		# transform labels to one-hot encoding
		self._classes = np.array(list(set(y)))

		n_classes = len(self._classes)
		y = self._onehot_encode(y)

		# shuffle training data
		x, y = self._shuffle(x, y)

		# initialize model
		self._initialize(n_input, n_classes)

		# initialize checkpoint saver
		saver = tf.train.Saver()

		# restore checkpoint if specified
		if self.load:
			saver.restore(self._session, "../checkpoints/dataset_nn")

		# perform training
		for epoch in range(self.epochs):
			# determine number of batches
			n_batches = int(n_samples / self.batch_size)
			avg_loss = 0

			# train on each batch
			for i in range(n_batches):
				# extract batch
				batch_x, batch_y = self._next_batch(x, y, self.batch_size, i)

				# compute loss
				_, loss, _ = self._session.run([self._optimizer, self._loss, self._model], feed_dict={ self._x: batch_x, self._y: batch_y })

				avg_loss += loss / n_batches

			# print results
			if self.verbose:
				print("epoch: %04d, lr: %0.6f, loss: %0.6f" % (epoch + 1, self._lr.eval(feed_dict=None, session=self._session), avg_loss))

		# save checkpoint if specified
		if self.save:
			saver.save(self._session, "../checkpoints/dataset_nn")



	def predict(self, x):
		return self._model.eval({ self._x: x }, session=self._session)



	def score(self, x, y):
		# transform labels to one-hot encoding
		y = self._onehot_encode(y)

		# compute accuracy of model on test data
		accuracy = tf.equal(tf.argmax(self._model, 1), tf.argmax(self._y, 1))
		accuracy = tf.reduce_mean(tf.cast(accuracy, "float"))

		return accuracy.eval({ self._x: x, self._y: y }, session=self._session)
