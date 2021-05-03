import tensorflow as tf
from tensorflow.keras import layers, metrics
from tensorflow.keras.layers.experimental import preprocessing as prepro
from tensorflow import feature_column
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import time
from sklearn.model_selection import train_test_split
import sklearn.preprocessing as skpp
import sklearn.utils


def nn1_quad(
    learning_rate=1e-3,
    weight_decay=1e-1,
    batch_size=256,
    epochs=3000,
    polygon_size=6,
):
    with tf.device('/cpu:0'):

        tmp_polygons = pd.read_csv("data/6-gon-quad-mesh-dataset.csv")
        polygons = tmp_polygons.copy()
        polygons.drop(polygons.columns[14:], axis=1, inplace=True)

        polygons = sklearn.utils.shuffle(polygons)
        # Split dataset into 70/15/15 training/test
        polygon_train = polygons.sample(frac=0.85, random_state=0)
        polygon_test = polygons.drop(polygon_train.index)

        train_features = polygon_train.copy()
        train_labels = train_features.pop('internal_node_count')

        test_features = polygon_test.copy()
        test_labels = test_features.pop('internal_node_count')

        polygon_model = tf.keras.Sequential([
            tf.keras.Input(shape=(13,)),
            layers.BatchNormalization(),

            layers.Dense(24),
            layers.Activation('relu'),
            layers.BatchNormalization(),
            layers.Dropout(0.5),

            layers.Dense(24),
            layers.Activation('relu'),
            layers.BatchNormalization(),
            layers.Dropout(0.5),

            layers.Dense(1),
        ])

        polygon_model.summary()
        polygon_model.compile(loss=tf.losses.MeanAbsoluteError(),
                              optimizer=tf.optimizers.Adam(
                                  learning_rate=learning_rate, decay=weight_decay
        ),
        )

        history = polygon_model.fit(train_features,
                                    train_labels,
                                    epochs=epochs,
                                    batch_size=batch_size,
                                    validation_split=0.18,
                                    verbose=2,
                                    )

        # Evaluate the model
        train_acc = polygon_model.evaluate(
            train_features, train_labels, verbose=0)
        test_acc = polygon_model.evaluate(
            test_features, test_labels, verbose=0)
        print('Training data loss: %.3f, Test data loss: %.3f' %
              (train_acc, test_acc))

        # plot history
        plt.plot(history.history['loss'], label='training loss')
        plt.plot(history.history['val_loss'], label='validation loss')
        plt.legend()
        plt.show()


if __name__ == '__main__':
    start_time = time.time()
    nn1_quad()
    print(f"Time elapsed: {round(time.time() - start_time, 5)} seconds")
