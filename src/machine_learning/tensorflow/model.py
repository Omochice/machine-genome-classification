import tensorflow as tf
from tensorflow.keras.layers import Input, Conv2D, Concatenate, Dropout, Dense
from tensorflow.keras.applications import MobileNetV2
from tensorflow.keras.utils import plot_model


class MyModel(tf.keras.Model):
    def __init__(self):
        super(MyModel, self).__init__()
        self.input1 = tf.keras.layers.Input()
        self.change_channel = tf.keras.layers.Conv2D(filter=3, kernel_size=3,
                                                     padding="same", activation="relu")
        self.mobilent = tf.keras.applications.MobileNetV2(include_top=False,
                                                          input_shape=(
                                                              128, 128, 3),
                                                          weights="imagenet",
                                                          pooling="max")


def construct_model(n_class: int):
    input1 = Input(shape=(128, 128, 1))
    change_channel = Conv2D(filters=3, kernel_size=3, padding="same",
                            activation="relu")(input1)
    mobilenet = MobileNetV2(include_top=False,
                            input_shape=(128, 128, 3),
                            weights="imagenet",
                            pooling="max")(change_channel)
    input2 = Input(shape=(5, ))
    concate = Concatenate()([mobilenet, input2])
    dense = Dense(128, activation="relu")(concate)
    dropout = Dropout(0.3)(dense)
    pred = Dense(n_class, activation="softmax")(dropout)
    return tf.keras.Model(inputs=[input1, input2], outputs=pred)


def show_model(model, dst: str = "model.png") -> None:
    model.summary()
    plot_model(model, to_file=dst)
    print(f"model img saved in {dst}")


if __name__ == "__main__":
    pass
