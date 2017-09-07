#/usr/bin/python

# multilayer perceptron neural network with softmax layer to classify
# GTEx data (30 classes)

import numpy as np 
import tensorflow as tf 
import setup_gtex

# Parameters
learning_rate = 0.001
training_epochs = 100
batch_size = 16
display_step = 1
beta = 0.01

# Network Parameters
n_hidden_1 = 8192 # 1st layer number of features
n_hidden_2 = 4096 # 2nd layer number of features
n_hidden_3 = 512 # 3rd layer number of features
#n_hidden_4 = 128 # 4th layer num neurons
n_input = 56238 # GTEx data input size
n_classes = 30 # GTEx total classes

# tf Graph input
x = tf.placeholder("float", [None, n_input])
y = tf.placeholder("float", [None, n_classes])

# Create model
def multilayer_perceptron(x, weights, biases):
    # Hidden layer with RELU activation
    #x = tf.nn.l1_normalize(x, 0)
    layer_1 = tf.add(tf.matmul(x, weights['h1']), biases['b1'])
    layer_1 = tf.nn.relu(layer_1)
    #layer_1 = tf.nn.dropout(layer_1, 0.50)
    # Hidden layer with RELU activation
    layer_2 = tf.add(tf.matmul(layer_1, weights['h2']), biases['b2'])
    layer_2 = tf.nn.relu(layer_2)
    #layer_2 = tf.nn.dropout(layer_2, 0.50)
    # Hidden layer with RELU activation
    layer_3 = tf.add(tf.matmul(layer_2, weights['h3']), biases['b3'])
    layer_3 = tf.nn.relu(layer_3)
    #layer_3 = tf.nn.dropout(layer_3, 0.50)
    #layer_3 = tf.nn.l2_normalize(layer_3, [0,1])


    # Hidden layer with Sigmoid activation
    #layer_4 = tf.add(tf.matmul(layer_3, weights['h4']), biases['b4'])
    #layer_4 = tf.nn.sigmoid(layer_4)
    # Output layer with linear activation
    out_layer = tf.matmul(layer_3, weights['out']) + biases['out']
    return out_layer

# Store layers weight & bias
weights = {
    'h1': tf.Variable(tf.random_normal([n_input, n_hidden_1])),
    'h2': tf.Variable(tf.random_normal([n_hidden_1, n_hidden_2])),
    'h3': tf.Variable(tf.random_normal([n_hidden_2, n_hidden_3])),
    #'h4': tf.Variable(tf.random_normal([n_hidden_3, n_hidden_4])),
    'out': tf.Variable(tf.random_normal([n_hidden_3, n_classes]))
}
biases = {
    'b1': tf.Variable(tf.random_normal([n_hidden_1])),
    'b2': tf.Variable(tf.random_normal([n_hidden_2])),
    'b3': tf.Variable(tf.random_normal([n_hidden_3])),
    #'b4': tf.Variable(tf.random_normal([n_hidden_4])),
    'out': tf.Variable(tf.random_normal([n_classes]))
}

# gather data
gtex = setup_gtex.GTEx('./datasets/GTEx_Data_30', './train_data', './test_data')

# Construct model
pred = multilayer_perceptron(x, weights, biases)

# Define loss and optimizer
result = tf.nn.softmax(pred)
cost = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(logits=pred, labels=y))
#regularizer = tf.nn.l2_loss(weights['h1'])  + tf.nn.l2_loss(weights['h2']) + tf.nn.l2_loss(weights['h3']) + tf.nn.l2_loss(weights['out']) \
#                + tf.nn.l2_loss(biases['b1']) + tf.nn.l2_loss(biases['b2']) + tf.nn.l2_loss(biases['b3']) + tf.nn.l2_loss(biases['out'])
#cost = tf.reduce_mean(cost)# + beta * regularizer)

optimizer = tf.train.AdamOptimizer(learning_rate=learning_rate).minimize(cost)
saver = tf.train.Saver()

# Initializing the variables
init = tf.global_variables_initializer()

# Launch the graph
sess = tf.Session()
sess.run(init)

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
        
    if epoch % display_step == 0:
        print("Epoch:", '%04d' % (epoch+1), "cost=", \
            "{:.9f}".format(avg_cost))
print("Optimization Finished!")
saver.save(sess, "./checkpoints/gtex_nn")

# Test model
correct_prediction = tf.equal(tf.argmax(pred, 1), tf.argmax(y, 1))
# Calculate accuracy
accuracy = tf.reduce_mean(tf.cast(correct_prediction, "float"))
print("Accuracy:", accuracy.eval({x: gtex.test.data, y: gtex.test.labels}, session=sess))


# res = tf.argmax(pred, 1)

# a =  sess.run(res, feed_dict={x: gtex.test.data})
