import tensorflow as tf
from tensorflow.keras.layers import Input, Conv2D, Concatenate, Dropout, Dense
from tensorflow.keras.applications import MobileNetV2, MobileNet
from tensorflow.keras.utils import plot_model


class MyModel(tf.keras.Model):
    def __init__(self):
        super(MyModel, self).__init__()
        self.input1 = tf.keras.layers.Input()
        self.change_channel = tf.keras.layers.Conv2D(filter=3,
                                                     kernel_size=3,
                                                     padding="same",
                                                     activation="relu")
        self.mobilent = tf.keras.applications.MobileNetV2(include_top=False,
                                                          input_shape=(128, 128, 3),
                                                          weights="imagenet",
                                                          pooling="max")


# Optunaでパラメータチューニングしたい
def construct_model(n_class: int) -> tf.keras.Model:
    input1 = Input(shape=(192, 192, 1))
    change_channel = Conv2D(filters=3, kernel_size=3, padding="same",
                            activation="relu")(input1)
    mobilenet = MobileNetV2(include_top=False,
                            input_shape=(192, 192, 3),
                            weights="imagenet",
                            pooling="max")(change_channel)
    input2 = Input(shape=(1, ))
    concate = Concatenate()([mobilenet, input2])
    dense1 = Dense(512, activation="relu")(concate)
    dropout1 = Dropout(0.3)(dense1)
    dense2 = Dense(512, activation="relu")(dropout1)
    dropout2 = Dropout(0.3)(dense2)
    pred = Dense(n_class, activation="softmax")(dropout2)
    return tf.keras.Model(inputs=[input1, input2], outputs=pred)


def show_model(model, dst: str = "model.png") -> None:
    model.summary()
    plot_model(model, to_file=dst, show_shapes=True)
    print(f"model img saved in {dst}")


if __name__ == "__main__":
    pass
