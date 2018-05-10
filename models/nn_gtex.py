#/usr/bin/python

# multilayer perceptron neural network with softmax layer to classify genetic data
import numpy as np 
import tensorflow as tf 
from sklearn import preprocessing
import sys, argparse
import os

import matplotlib.pyplot as plt

def confusion_heatmap(conf_arr):
    norm_conf = preprocessing.normalize(conf_arr, axis=1, norm='l1')

    fig = plt.figure(figsize=(17,10))
    plt.clf()
    ax = fig.add_subplot(111)
    ax.set_aspect(1)

    res = ax.imshow(norm_conf, cmap=plt.cm.jet, 
            interpolation='nearest')

    width, height = conf_arr.shape

    # for x in xrange(width):
    #     for y in xrange(height):
    #         ax.annotate(str(conf_arr[x][y]), xy=(y, x), 
    #                     horizontalalignment='center',
    #                     verticalalignment='center')

    cb = fig.colorbar(res)
    plt.xticks(range(width), np.arange(0,conf_arr.shape[0]), rotation='vertical')
    plt.yticks(range(height), np.arange(0,conf_arr.shape[0]))
    plt.show()
    # plt.savefig('confusion_matrix.png', format='png')


class MLP:
    def __init__(self, lr=0.001, epochs=75, n_layers=3, h_units=[512,512,512], \
        act_funcs=["relu", "relu", "relu"], batch_size=16, disp_step=1, n_input=56238, \
        n_classes=53, dropout=0, load=0, save = 0, confusion=0, verbose=0):
        
        self.lr = lr
        self.epochs = epochs
        self.n_layers = n_layers
        self.h_units = h_units
        self.act_funcs = act_funcs
        self.batch_size = batch_size
        self.display_step = disp_step
        self.n_input = n_input
        self.n_classes = n_classes
        self.load = load
        self.save = save
        self.dropout = dropout
        self.confusion = confusion
        self.verbose = verbose

        os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' 

    # Create model
    def multilayer_perceptron(self, x, weights, biases):
        
        layer = x
        for i in xrange(1, self.n_layers + 1):
            w = 'h' + str(i)
            b = 'b' + str(i)
            
            layer = tf.add(tf.matmul(layer, weights[w]), biases[b])
            
            if self.act_funcs[i - 1] == "relu":
                layer = tf.nn.relu(layer)
            elif self.act_funcs[i - 1] == "sigmoid":
                layer = tf.nn.sigmoid(layer)

            if self.dropout:
                layer = tf.nn.dropout(layer, 0.75)

        out_layer = tf.add(tf.matmul(layer, weights['out']), biases['out'])

        return out_layer


    def run(self, gtex):

        tf.reset_default_graph()

        x = tf.placeholder("float", [None, self.n_input])
        y = tf.placeholder("float", [None, self.n_classes])

        units = [self.n_input]
        for i in self.h_units:
            units.append(i)
        units.append(self.n_classes)

        weights = {}
        biases = {}
        for i in xrange(1, self.n_layers + 1):
            w = 'h' + str(i)
            b = 'b' + str(i)
            weights[w] = tf.get_variable(w, shape=[units[i - 1], units[i]], initializer=tf.contrib.layers.xavier_initializer())
            biases[b] = tf.get_variable(b, shape=[units[i]], initializer=tf.contrib.layers.xavier_initializer())
        
        weights['out'] = tf.get_variable('out_w', shape=[self.h_units[-1], self.n_classes], initializer=tf.contrib.layers.xavier_initializer())
        biases['out'] = tf.get_variable('out_b', shape=[self.n_classes], initializer=tf.contrib.layers.xavier_initializer())

        # preprocess data
        maxabsscaler = preprocessing.MaxAbsScaler()
        gtex.train.data = maxabsscaler.fit_transform(gtex.train.data)
        gtex.test.data = maxabsscaler.fit_transform(gtex.test.data)
        
        # Construct model
        pred = self.multilayer_perceptron(x, weights, biases)

        # Define loss and optimizer
        g_step = tf.Variable(0, trainable=False)

        # set learning rate that will be decayed
        starter_learning_rate = self.lr

        learning_rate = tf.train.exponential_decay(starter_learning_rate, global_step=g_step, decay_steps=500, decay_rate=0.96, staircase=True)

        result = tf.nn.softmax(pred)
        cost = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(logits=pred, labels=y))

        optimizer = tf.train.AdamOptimizer(learning_rate=learning_rate).minimize(cost, global_step=g_step)
        saver = tf.train.Saver()

        # Initializing the variables
        init = tf.global_variables_initializer()

        # Launch the graph
        sess = tf.Session()
        sess.run(init)

        if self.load:
            saver.restore(sess, '../checkpoints/gtex_nn')

        # Training cycle
        for epoch in range(self.epochs):
            avg_cost = 0.
            total_batch = int(gtex.train.num_examples/self.batch_size)
            #gtex.shuffle()
            # Loop over all batches
            for i in range(total_batch):
                batch_x, batch_y = gtex.train.next_batch(self.batch_size, i)

                # Run optimization op (backprop) and cost op (to get loss value)
                _, c, r = sess.run([optimizer, cost, result], feed_dict={x: batch_x,
                                                              y: batch_y})
                # Compute average loss
                avg_cost += c / total_batch

            if self.verbose:
                print("Epoch:", '%04d' % (epoch+1), "Learning Rate: ", '%5f' % (learning_rate.eval(feed_dict=None, session=sess)), "cost=", "{:.9f}".format(avg_cost))

        if self.save:
            saver.save(sess, "../checkpoints/gtex_nn")

        # Test model
        correct_prediction = tf.equal(tf.argmax(pred, 1), tf.argmax(y, 1))
        
        # Calculate accuracy
        accuracy = tf.reduce_mean(tf.cast(correct_prediction, "float"))

        if self.confusion:
            # generate confusion matrices for brain data and gtex data
            temp = pred.eval({x: gtex.test.data}, session=sess)
            preds = np.argmax(temp, 1)
            labs = np.argmax(gtex.test.labels, 1)
            cm = tf.confusion_matrix(labs, preds, num_classes=self.n_classes)
            mycm = cm.eval(feed_dict=None, session=sess)
            print mycm
            np.savetxt('./confusion_matrix_gtex.txt', mycm, fmt='%4d', delimiter=' ')

        # calculate accuracy that will be returned
        acc = accuracy.eval({x: gtex.test.data, y: gtex.test.labels}, session=sess)

        sess.close()

        return acc
