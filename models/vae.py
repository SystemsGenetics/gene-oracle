#/usr/bin/python

'''
	This script contains the code to train and evaluate a variational autoencoder on
	genetic RNA data (but can be applicable to other datatypes)

	Protypes:

	Todo:

'''

import tensorflow as tf
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os


markers = {',': 'pixel', 'o': 'circle','*': 'star', 'v': 'triangle_down',
           '^': 'triangle_up', '<': 'triangle_left', '>': 'triangle_right', 
           '1': 'tri_down', '2': 'tri_up', '3': 'tri_left', '4': 'tri_right', 
           '8': 'octagon', 's': 'square', 'p': 'pentagon', 
           'h': 'hexagon1', 'H': 'hexagon2', '+': 'plus', 'x': 'x', '.': 'point', 
           'D': 'diamond', 'd': 'thin_diamond', '|': 'vline', '_': 'hline',
           'P': 'plus_filled', 'X': 'x_filled', 0: 'tickleft', 
           1: 'tickright', 2: 'tickup', 3: 'tickdown', 4: 'caretleft', 5: 'caretright',
           6: 'caretup', 7: 'caretdown', 8: 'caretleftbase', 9: 'caretrightbase', 10: 'caretupbase',
           11: 'caretdownbase', 'None': 'nothing', None: 'nothing', ' ': 'nothing', '': 'nothing'}
markers_keys = list(markers.keys())[:20]

font = {'family' : 'normal',
         'weight' : 'bold',
         'size'   : 30}

plt.rc('font', **font)

sns.set_style("ticks")

colors = ["windows blue", "amber", 
          "greyish", "faded green", 
          "dusty purple","royal blue","lilac",
          "salmon","bright turquoise",
          "dark maroon","light tan",
          "orange","orchid",
          "sandy","topaz",
          "fuchsia","yellow",
          "crimson","cream"
          ]
current_palette = sns.xkcd_palette(colors)


def print_2D( points,label):
    '''
    points: N_samples * 2
    label: (int) N_samples
    id_map: map label id to its name
    '''  
    fig = plt.figure()
    #current_palette = sns.color_palette("RdBu_r", max(label)+1)
    n_cell,_ = points.shape
    if n_cell > 500:
        s = 10
    else:
        s = 20

    id_map = {i:str(i) for i in range(10)}

    label = np.argmax(label, axis=1)
    
    ax = plt.subplot(111)
    print( np.unique(label) )
    for i in np.unique(label):
        ax.scatter( points[label==i,0], points[label==i,1], c=current_palette[i], label=id_map[i], s=s,marker=markers_keys[i] )
    box = ax.get_position()
    # ax.set_position([box.x0, box.y0 + box.height * 0.1,
    #                  box.width, box.height * 0.9])
        
    ax.legend(scatterpoints=1,loc='upper center',
              bbox_to_anchor=(0.5,-0.08),ncol=6,
              fancybox=True,
              prop={'size':8}
              )
    sns.despine()
    return fig

def plot(samples):
    fig = plt.figure(figsize=(4, 4))
    gs = gridspec.GridSpec(4, 4)
    gs.update(wspace=0.05, hspace=0.05)

    for i, sample in enumerate(samples):
        ax = plt.subplot(gs[i])
        plt.axis('off')
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_aspect('equal')
        plt.imshow(sample.reshape(28, 28), cmap='Greys_r')

    return fig

class VAE:
	def __init__(self, latent_dim=2, lr=0.001, epochs=75, h_units=[512,256,128], \
				batch_size=64, n_input=56238, dropout=0, load=0, save=0, verbose=0):
		self.latent_dim = latent_dim
		self.lr = lr
		self.epochs = epochs
		self.h_units = h_units
		self.batch_size = batch_size
		self.n_input = n_input
		self.load = load
		self.save = save
		self.dropout = dropout
		self.verbose = verbose

		os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'


	# USAGE:
	# 		- encoder network for vae
	# PARAMS:
	#	x: input data sample
	#	h_hidden: LIST of num. neurons per hidden layer
	def inference_network(self, x):
		with tf.variable_scope('encoder', reuse=tf.AUTO_REUSE):
			layer = x
			for n in self.h_units:
				layer = tf.layers.dense(inputs=layer, units=n, activation=tf.nn.relu)

			gaussian_params = tf.layers.dense(inputs=layer, units=self.latent_dim * 2, activation=None)

			mu = gaussian_params[:, :self.latent_dim]

			# std dev must be positive... parametrize with a softplus
			sigma = tf.nn.softplus(gaussian_params[:, self.latent_dim:])

			return mu, sigma


	# USAGE:
	# 		- decoder network for vae
	# PARAMS:
	#	z: input latent variable
	#	n_hidden: LIST of num. neurons per hidden layer
	def generative_network(self, z):
		with tf.variable_scope('decoder', reuse=tf.AUTO_REUSE):
			layer = z
			for n in self.h_units:
				layer = tf.layers.dense(inputs=layer, units=n, activation=tf.nn.relu)

			# use a sigmoid activation to get the data between 0 and 1
			logits = tf.layers.dense(inputs=layer, units=self.n_input, activation=None)
			probs = tf.nn.sigmoid(logits)

		return probs, logits



	# USAGE:
	# 		- perform the reparameterization trick for vae
	# PARAMS:
	#	mean: mean produced by inference network
	#	sigma: sigma produced by inference network
	def reparameterize(self, mean, logvar):
		eps = tf.random_normal(shape=tf.shape(mean))
		q_z = mean + tf.exp(logvar * 0.5) * eps	
		return q_z


	def train(self, dataset):
		# define placeholders for input data
		x = tf.placeholder("float", [None, self.n_input])
		z = tf.placeholder(tf.float32, shape=[None, self.latent_dim])

		# # preprocess data
		# maxabsscaler = preprocessing.MaxAbsScaler()
		# dataset.train.data = maxabsscaler.fit_transform(dataset.train.data)
		# dataset.test.data = maxabsscaler.fit_transform(dataset.test.data)

		# first run the inference model to produce the gaussian parameters
		q_mu, q_sigma = self.inference_network(x=x)

		q_z = self.reparameterize(q_mu, q_sigma)

		_, x_logit = self.generative_network(q_z)

		X_samples, _ = self.generative_network(z)

		# define losses
		recon_loss = tf.reduce_sum(tf.nn.sigmoid_cross_entropy_with_logits(logits=x_logit, labels=x), axis=1)
		kl = 0.5 * tf.reduce_sum(tf.exp(q_sigma) + tf.square(q_mu) - 1. - q_sigma, axis=1)

		ELBO = tf.reduce_mean(recon_loss + kl)

		loss = ELBO

		# optimizer
		optimizer = tf.train.AdamOptimizer(learning_rate=self.lr).minimize(loss)

		saver = tf.train.Saver()

		# Initializing the variables
		init = tf.global_variables_initializer()

		sess = tf.Session()
		sess.run(init)

		for epoch in range(1, self.epochs + 1):
			avg_cost = 0.
			total_batch = int(dataset.train.num_examples/self.batch_size)

			for i in range(total_batch):
				batch_x, batch_y = dataset.train.next_batch(self.batch_size, i)
				_, c, rcl, kll = sess.run([optimizer, loss, recon_loss, kl], feed_dict={x: batch_x})

				avg_cost += c / total_batch

			print("Epoch:", '%04d' % (epoch), "cost=", "{:.9f}".format(avg_cost))

		saver.save(sess, "/tmp/model.ckpt")
		sess.close() 


	# run inference on the inference network Q(z | X)
	def infer(self, data):
		x = tf.placeholder("float", [None, self.n_input])
		q_mu, q_sigma = self.inference_network(x=x)

		saver = tf.train.Saver()

		sess = tf.Session()
		saver.restore(sess, "/tmp/model.ckpt")

		mu, sig = sess.run([q_mu, q_sigma], feed_dict={x: data})

		sess.close()

		return mu, sig


	# run inference on the generative network P(X | z)
	def generate(self):
		z = tf.placeholder(tf.float32, shape=[None, self.latent_dim])

		X_samples, _ = self.generative_network(z)

		saver = tf.train.Saver()

		sess = tf.Session()

		saver.restore(sess, "/tmp/model.ckpt")

		#mu, sig = 
		xsamp = sess.run(X_samples, feed_dict={z: np.random.randn(16, self.latent_dim)})

		fig = plot(xsamp)
		plt.savefig('out/{}.png'.format(str(0).zfill(3)), bbox_inches='tight')
		plt.close(fig)

		sess.close()











