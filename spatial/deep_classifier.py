import tensorflow as tf
from tensorflow.keras import optimizers, regularizers
from tensorflow.keras.models import Model, load_model, Sequential
from tensorflow.keras.layers import Dense, Activation, BatchNormalization, Dropout
import numpy as np


class CLASSIFIER:
    def __init__(self, input_size, class_num=2, path=''):
        self.input_size = input_size
        self.classifier = None
        self.initializers = "glorot_uniform"
        self.optimizer = optimizers.Adam(lr=0.01)
        self.validation_split = 0.1
        self.class_num = class_num
        self.dropout_rate = 0.05

    def build(self):
        model = Sequential()
        model.add(Dense(64, activation='relu', input_shape=(self.input_size,)))
        model.add(Dropout(rate=self.dropout_rate))
        model.add(Dense(32, activation='relu'))
        model.add(Dropout(rate=self.dropout_rate))
        model.add(Dense(self.class_num, activation='softmax'))
        self.classifier = model

    def compile(self):
        self.classifier.compile(optimizer=self.optimizer, loss='categorical_crossentropy', metrics=['accuracy'])
        self.classifier.summary()

    def train(self, x, label, batch_size=100, epochs=300):
        history = self.classifier.fit(x, label,
                                    epochs=epochs, batch_size=batch_size,
                                    validation_split=self.validation_split, shuffle=True)
        return history

    def prediction(self, x):
        label = self.classifier.predict(x)
        label = np.argmax(label, axis=1)
        return label
