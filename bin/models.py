import numpy as np
import sklearn.base
import sys

from tensorflow import keras



class TensorflowMLP(sklearn.base.BaseEstimator):

    def __init__(self,
                 hidden_layer_sizes=[512, 256, 128],
                 activation='relu',
                 lr=0.001,
                 epochs=50,
                 batch_size=32,
                 verbose=False):

        # save attributes
        self.hidden_layer_sizes = hidden_layer_sizes
        self.activation = activation
        self.lr = lr
        self.epochs = epochs
        self.batch_size = batch_size
        self.verbose = verbose

    def _build(self, n_inputs, n_classes):
        # save attributes
        self.n_inputs = n_inputs
        self.n_classes = n_classes

        # delete previous model
        keras.backend.clear_session()

        # define input layer
        x_input = keras.Input(shape=n_inputs)

        # define hidden layers
        x = x_input
        for units in self.hidden_layer_sizes:
            x = keras.layers.Dense(units=units, activation=self.activation)(x)

        # define output layer
        y_output = keras.layers.Dense(units=n_classes, activation='softmax')(x)

        # define model
        model = keras.models.Model(x_input, y_output)

        # define training parameters
        model.compile(
            optimizer=keras.optimizers.Adam(learning_rate=self.lr),
            loss='categorical_crossentropy',
            metrics=['accuracy'])

        # save model
        self._model = model

    def _onehot_encode(self, y):
        return keras.utils.to_categorical(y, num_classes=self.n_classes)

    def fit(self, x, y):
        # build model
        n_inputs = x.shape[1]
        n_classes = len(set(y))

        self._build(n_inputs, n_classes)

        # transform labels to one-hot encoding
        y = self._onehot_encode(y)

        # train model
        return self._model.fit(
            x=x,
            y=y,
            batch_size=self.batch_size,
            epochs=self.epochs,
            verbose=self.verbose)

    def predict(self, x):
        # get predictions
        y_pred = self._model.predict(x)

        # decode predictions from one-hot
        y_pred = np.argmax(y_pred, axis=1)

        return y_pred

    def score(self, x, y):
        # transform labels to one-hot encoding
        y = self._onehot_encode(y)

        # evaluate model
        loss, accuracy = self._model.evaluate(x, y)

        # return accuracy
        return accuracy
