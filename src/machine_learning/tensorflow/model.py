import numpy as np
from tensorflow.keras.utils import plot_model
from pathlib import Path
import tensorflow as tf


class MyModel(tf.keras.Model):
    def __init__(self) -> None:
        super(MyModel, self).__init__()
        self.flatten = tf.keras.layers.Flatten()
        self.fc1 = tf.keras.layers.Dense(128, activation="relu")
        self.fc2 = tf.keras.layers.Dense(128, activation="softmax")
        self.dropout = tf.keras.layers.Dropout(0.2)

    def call(self, x, training=False):
        x = self.flatten(x)
        x = self.fc1(x)
        x = self.dropout(x, training=training)
        x = self.fc2(x)
        return x

    def build_graph(self, dim):
        x = tf.keras.layers.Input(shape=(dim))
        return tf.keras.Model(inputs=[x], outputs=self.call(x))


if __name__ == "__main__":
    model = MyModel()
    model.compile(optimizer="sgd",
                  loss="categorical_crossentropy",
                  metrics=["accuracy"])
    model.build((None, 20, 20))
    model.summary()
    plot_model(model.build_graph((None, 20, 20)), to_file="model.png")
