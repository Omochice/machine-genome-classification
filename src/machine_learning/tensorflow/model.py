import numpy as np
from tensorflow.keras.utils import plot_model
from pathlib import Path
import tensorflow as tf
from tensorflow.keras.applications import MobileNet


class MyModel(tf.keras.Model):
    def __init__(self, img_dim: tuple, f_dim) -> None:
        super(MyModel, self).__init__()
        self.image_input = tf.keras.layers.Input(shape=(img_dim))
        self.linear_input = tf.keras.layers.Input(shape=(100, ))
        self.concate = tf.keras.layers.Concatenate()
        self.mobilenet = MobileNet(include_top=False,
                                   input_shape=(128, 128, 3),
                                   weights="imagenet",
                                   pooling="max")
        self.change_chanel = tf.keras.layers.Conv2D(filters=3,
                                                    kernel_size=3,
                                                    padding="same",
                                                    activation="relu")
        self.flatten = tf.keras.layers.Flatten()
        self.fc1 = tf.keras.layers.Dense(512, activation="relu")
        self.fc2 = tf.keras.layers.Dense(128, activation="softmax")
        self.dropout = tf.keras.layers.Dropout(0.2)

    def call(self, inputs, training=False):
        x = self.change_chanel(inputs[0])
        x = self.mobilenet(x)
        x = self.concate([x, inputs[1]])
        x = self.flatten(x)
        x = self.fc1(x)
        x = self.dropout(x, training=training)
        x = self.fc2(x)
        return x

    def build_graph(self, dims):
        x = tf.keras.layers.Input(shape=(dims[0]))
        y = tf.keras.layers.Input(shape=(dims[1]))
        return tf.keras.Model(inputs=[x, y], outputs=self.call((x, y)))


if __name__ == "__main__":
    shape1 = (None, 128, 128, 1)
    shape2 = (None, 2)
    model = MyModel(shape1, shape2)
    model.compile(optimizer="sgd",
                  loss="categorical_crossentropy",
                  metrics=["accuracy"])
    model.build(input_shape=[shape1, shape2])
    model.summary()
    plot_model(model.build_graph((shape1, shape2)), to_file="model.png")
