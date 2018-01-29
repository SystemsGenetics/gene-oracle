#/usr/bin/python

# multilayer perceptron neural network with softmax layer to classify
# GTEx data (30 classes)

import sys
sys.path.insert(0, '../data_scripts')

import numpy as np 
import tensorflow as tf 
from sklearn import preprocessing
import setup_gtex
import sys, argparse
import os

class MLP:
    def __init__(self, lr=0.001, epochs=100, n_h1=512, n_h2=512, n_h3=512, batch_size=16, \
        disp_step=1, n_input=56238, n_classes=53, beta=0.01, load=0, confusion=0, verbose=0):
        
        self.lr = lr
        self.epochs = epochs
        self.n_hidden_1 = n_h1
        self.n_hidden_2 = n_h2
        self.n_hidden_3 = n_h3
        self.batch_size = batch_size
        self.display_step = disp_step
        self.n_input = n_input
        self.n_classes = n_classes
        self.beta = beta
        self.load = load
        self.confusion = confusion
        self.verbose = verbose

        os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' 

    # Create model
    def multilayer_perceptron(self, x, weights, biases):
        # Hidden layer with RELU activation
        #x = tf.nn.l1_normalize(x, 0)
        layer_1 = tf.add(tf.matmul(x, weights['h1']), biases['b1'])
        layer_1 = tf.nn.relu(layer_1)
        #layer_1 = tf.nn.dropout(layer_1, 0.75)
        # Hidden layer with RELU activation
        layer_2 = tf.add(tf.matmul(layer_1, weights['h2']), biases['b2'])
        layer_2 = tf.nn.relu(layer_2)
        #layer_2 = tf.nn.dropout(layer_2, 0.75)
        # Hidden layer with RELU activation
        layer_3 = tf.add(tf.matmul(layer_2, weights['h3']), biases['b3'])
        layer_3 = tf.nn.sigmoid(layer_3)
        #layer_3 = tf.nn.dropout(layer_3, 0.75)
        #layer_3 = tf.nn.l2_normalize(layer_3, [0,1])


        # Hidden layer with Sigmoid activation
        #layer_4 = tf.add(tf.matmul(layer_3, weights['h4']), biases['b4'])
        #layer_4 = tf.nn.sigmoid(layer_4)
        # Output layer with linear activation
        out_layer = tf.matmul(layer_3, weights['out']) + biases['out']
        return out_layer


    def run(self, gtex):

        tf.reset_default_graph()

        x = tf.placeholder("float", [None, self.n_input])
        y = tf.placeholder("float", [None, self.n_classes])

        # Store layers weight & bias
        weights = {
            'h1': tf.get_variable("h1", shape=[self.n_input, self.n_hidden_1], initializer=tf.contrib.layers.xavier_initializer()),
            'h2': tf.get_variable("h2", shape=[self.n_hidden_1, self.n_hidden_2], initializer=tf.contrib.layers.xavier_initializer()),
            'h3': tf.get_variable("h3", shape=[self.n_hidden_2, self.n_hidden_3], initializer=tf.contrib.layers.xavier_initializer()),
            #'h4': tf.Variable(tf.random_normal([n_hidden_3, n_hidden_4])),
            'out': tf.get_variable("out_w", shape=[self.n_hidden_3, self.n_classes], initializer=tf.contrib.layers.xavier_initializer())
        }
        biases = {
            'b1': tf.get_variable("b1", shape=[self.n_hidden_1], initializer=tf.contrib.layers.xavier_initializer()),
            'b2': tf.get_variable("b2", shape=[self.n_hidden_2], initializer=tf.contrib.layers.xavier_initializer()),
            'b3': tf.get_variable("b3", shape=[self.n_hidden_3], initializer=tf.contrib.layers.xavier_initializer()),
            #'b4': tf.Variable(tf.random_normal([n_hidden_4])),
            'out': tf.get_variable("out_b", shape=[self.n_classes], initializer=tf.contrib.layers.xavier_initializer())
        }

        # gather data
        #gtex = setup_gtex.GTEx('../datasets/GTEx_Data', '../train_data', '../test_data', self.n_input)

        # preprocess data
        maxabsscaler = preprocessing.MaxAbsScaler()
        gtex.train.data = maxabsscaler.fit_transform(gtex.train.data)
        gtex.test.data = maxabsscaler.fit_transform(gtex.test.data)
        

        # Construct model
        pred = self.multilayer_perceptron(x, weights, biases)

        # Define loss and optimizer
        g_step = tf.Variable(0, trainable=False)

        starter_learning_rate = self.lr

        learning_rate = tf.train.exponential_decay(starter_learning_rate, global_step=g_step, decay_steps=500, decay_rate=0.96, staircase=True)

        result = tf.nn.softmax(pred)
        cost = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(logits=pred, labels=y))
        # l1_regularizer = tf.contrib.layers.l1_regularizer(scale=0.005, scope=None)
        # w = tf.trainable_variables()
        # regularize_penalty = tf.contrib.layers.apply_regularization(l1_regularizer, w) 
        # cost = tf.reduce_mean(cost + regularize_penalty)

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
                
        #     if epoch % display_step == 0:
        #         print("Epoch:", '%04d' % (epoch+1), "Learning Rate: ", '%5f' % (learning_rate.eval(feed_dict=None, session=sess)), "cost=", "{:.9f}".format(avg_cost))
        # print("Optimization Finished!")
        #saver.save(sess, "../checkpoints/gtex_nn")

        # Test model
        correct_prediction = tf.equal(tf.argmax(pred, 1), tf.argmax(y, 1))
        
        # Calculate accuracy
        accuracy = tf.reduce_mean(tf.cast(correct_prediction, "float"))

        if self.confusion:
            # generate confusion matrices for brain data and gtex data
            temp = pred.eval({x: gtex.test.data}, session=sess)
            preds = np.argmax(temp, 1)
            labs = np.argmax(gtex.test.labels, 1)
            cm = tf.confusion_matrix(labs, preds, num_classes=args.n_classes)
            mycm = cm.eval(feed_dict=None, session=sess)
            np.savetxt('./confusion_matrix_gtex', mycm, fmt='%4d', delimiter=' ')

            temp = pred.eval({x: bd}, session=sess)
            preds = np.argmax(temp, 1)
            labs = np.argmax(labels, 1)
            cm = tf.confusion_matrix(labs, preds, num_classes=args.n_classes)
            mycm = cm.eval(feed_dict=None, session=sess)
            np.savetxt('./confusion_matrix_brain', mycm, fmt='%4d', delimiter=' ')

        acc = accuracy.eval({x: gtex.test.data, y: gtex.test.labels}, session=sess)

        sess.close()

        return acc