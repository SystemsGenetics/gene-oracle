#/usr/bin/python

# multilayer perceptron neural network with softmax layer to classify
# GTEx data (30 classes)

import numpy as np 
import tensorflow as tf 
from sklearn import preprocessing
import setup_gtex
import sys, argparse
import os

def main():
    parser = argparse.ArgumentParser(description='Neural network to classify genetic data')
    parser.add_argument('--lr', help='learning rate', type=float, default=0.0001)
    parser.add_argument('--epochs', help='no. of training epoch', type=int, default=100)
    parser.add_argument('--h1', help='no. of neurons in hidden layer 1', type=int, default=512)
    parser.add_argument('--h2', help='no. of neurons in hidden layer 2', type=int, default=512)
    parser.add_argument('--h3', help='no. of neurons in hidden layer 3', type=int, default=512)
    parser.add_argument('--batch_size', help='batch size', type=int, default=16)
    parser.add_argument('--display_step', help='print updates after this many steps', type=int, default=10)
    parser.add_argument('--n_input', help='number of input features', type=int, default=56238)
    parser.add_argument('--n_classes', help='number of classes', type=int, default=30)
    parser.add_argument('--beta', help='hyperparemeter for l1 regularization of weights', type=float, default=0.01)
    parser.add_argument('--load', help='load weights from previous run', type=bool, default=0)
    parser.add_argument('--confusion', help='generate confusion matrix (1) or no (0)', type=bool, default=0)
    args = parser.parse_args()

    # Parameters
    learning_rate = args.lr
    training_epochs = args.epochs
    batch_size = args.batch_size
    display_step = args.display_step
    beta = args.beta

    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

    # Network Parameters
    n_hidden_1 = args.h1 # 1st layer number of features
    n_hidden_2 = args.h2 # 2nd layer number of features
    n_hidden_3 = args.h3 # 3rd layer number of features
    n_hidden_4 = 256
    n_hidden_5 = 256
    n_hidden_6 = 256
    #n_hidden_4 = 256 # 4th layer num neurons
    n_input = args.n_input # GTEx data input size
    n_classes = args.n_classes # GTEx total classes

    # begin computational graph
    # tf Graph input
    x = tf.placeholder("float", [None, n_input])
    y = tf.placeholder("float", [None, n_classes])

    # Create model
    def multilayer_perceptron(x, weights, biases):
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
        layer_3 = tf.add(layer_1, layer_3)
        layer_3 = tf.nn.relu(layer_3)
        #layer_3 = tf.nn.dropout(layer_3, 0.75)
        #layer_3 = tf.nn.l2_normalize(layer_3, [0,1])

        # Hidden layer with Sigmoid activation
        layer_4 = tf.add(tf.matmul(layer_3, weights['h4']), biases['b4'])
        layer_4 = tf.nn.relu(layer_4)

        layer_5 = tf.add(tf.matmul(layer_4, weights['h5']), biases['b5'])
        layer_5 = tf.nn.relu(layer_5)

        layer_6 = tf.add(tf.matmul(layer_5, weights['h6']), biases['b6'])
        layer_6 = tf.add(layer_4, layer_6)
        layer_6 = tf.nn.relu(layer_6)

        # Output layer with linear activation
        out_layer = tf.matmul(layer_6, weights['out']) + biases['out']
        return out_layer

    # Store layers weight & bias
    weights = {
        'h1': tf.Variable(tf.random_normal([n_input, n_hidden_1])),
        'h2': tf.Variable(tf.random_normal([n_hidden_1, n_hidden_2])),
        'h3': tf.Variable(tf.random_normal([n_hidden_2, n_hidden_3])),
        'h4': tf.Variable(tf.random_normal([n_hidden_3, n_hidden_4])),
        'h5': tf.Variable(tf.random_normal([n_hidden_4, n_hidden_5])),
        'h6': tf.Variable(tf.random_normal([n_hidden_5, n_hidden_6])),
        #'h4': tf.Variable(tf.random_normal([n_hidden_3, n_hidden_4])),
        'out': tf.Variable(tf.random_normal([n_hidden_6, n_classes]))
    }
    biases = {
        'b1': tf.Variable(tf.random_normal([n_hidden_1])),
        'b2': tf.Variable(tf.random_normal([n_hidden_2])),
        'b3': tf.Variable(tf.random_normal([n_hidden_3])),
        'b4': tf.Variable(tf.random_normal([n_hidden_4])),
        'b5': tf.Variable(tf.random_normal([n_hidden_5])),
        'b6': tf.Variable(tf.random_normal([n_hidden_6])),
        #'b4': tf.Variable(tf.random_normal([n_hidden_4])),
        'out': tf.Variable(tf.random_normal([n_classes]))
    }

    # gather data
    #print('loading gtex data...')
    gtex = setup_gtex.GTEx('./datasets/GTEx_Data', './train_data', './test_data', args.n_input)

    #print('hidden layers: ' + str(args.h1) + 'x' + str(args.h2) + 'x' + str(args.h3))
    #print('epochs:        ' + str(args.epochs))
    #print('learning rate: ' + str(args.lr)) 
    #print('load:          ' + str(args.load))

    # preprocess data
    maxabsscaler = preprocessing.MaxAbsScaler()
    gtex.train.data = maxabsscaler.fit_transform(gtex.train.data)
    gtex.test.data = maxabsscaler.fit_transform(gtex.test.data)
    #gtex.train.data = preprocessing.normalize(gtex.train.data)
    #gtex.test.data = preprocessing.normalize(gtex.test.data)

    # Construct model
    pred = multilayer_perceptron(x, weights, biases)

    g_step = tf.Variable(0, trainable=False)
    starter_learning_rate = args.lr
    learning_rate = tf.train.exponential_decay(starter_learning_rate, global_step=g_step, decay_steps=500, decay_rate=0.96, staircase=True)

    # Define loss and optimizer
    result = tf.nn.softmax(pred)
    cost = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(logits=pred, labels=y))
    l1_regularizer = tf.contrib.layers.l1_regularizer(scale=0.005, scope=None)
    w = tf.trainable_variables()
    regularize_penalty = tf.contrib.layers.apply_regularization(l1_regularizer, w) 
    cost = tf.reduce_mean(cost + regularize_penalty)

    optimizer = tf.train.AdamOptimizer(learning_rate=learning_rate).minimize(cost, globale_step=g_step)
    saver = tf.train.Saver()

    # Initializing the variables
    init = tf.global_variables_initializer()

    # Launch the graph
    sess = tf.Session()
    sess.run(init)

    if args.
        saver.restore(sess, './checkpoints/gtex_nn')

    # Training cycle
    for epoch in range(training_epochs):
        avg_cost = 0.
        total_batch = int(gtex.train.num_examples/batch_size)
        # Loop over all batches
        for i in range(total_batch):
            batch_x, batch_y = gtex.train.next_batch(batch_size, i)
            # Run optimization op (backprop) and cost op (to get loss value)
            _, c, r = sess.run([optimizer, cost, result], feed_dict={x: batch_x,
                                                          y: batch_y})
            # Compute average loss
            avg_cost += c / total_batch
            
        #if epoch % display_step == 0:
            #print("Epoch:", '%04d' % (epoch+1), "cost=", \
                #"{:.9f}".format(avg_cost))
    #print("Optimization Finished!")
    saver.save(sess, "./checkpoints/gtex_nn")

    # Test model
    correct_prediction = tf.equal(tf.argmax(pred, 1), tf.argmax(y, 1))
    # Calculate accuracy
    accuracy = tf.reduce_mean(tf.cast(correct_prediction, "float"))
    print(str(accuracy.eval({x: gtex.test.data, y: gtex.test.labels}, session=sess)))
    #bd = np.transpose(np.load('./datasets/braincell_data_flt32.npy'))
    #bd = maxabsscaler.fit_transform(bd)
    #labels = np.zeros((331,30))
    #labels[:,5] = 1
    #print("Braincell Accuracy: ", accuracy.eval({x:bd, y:labels}, session=sess))
    #print(result.eval({x: gtex.test.data[0:3]}, session=sess))

    if args.confusion:
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

    # res = tf.argmax(pred, 1)

    # a =  sess.run(res, feed_dict={x: gtex.test.data})

if __name__ == "__main__":
    main()
