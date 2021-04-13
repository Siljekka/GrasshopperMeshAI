import tensorflow as tf
from tensorflow.keras import layers, metrics
from tensorflow import feature_column
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import time
from sklearn.model_selection import train_test_split
import sklearn.preprocessing as skpp
import sklearn.utils

# print("Num GPUs Available: ", len(tf.config.list_physical_devices('GPU')))

# import data from csv, https://www.tensorflow.org/guide/data#consuming_csv_data
def nn1():
    learning_rate = 1e-4
    weight_decay = 1e-1
    batch_size = 512
    epochs = 3000
    polygon_size = 6

    with tf.device('/cpu:0'):

        polygons = pd.read_csv("data/6-gon-mesh-dataset.csv")

        polygons = sklearn.utils.shuffle(polygons)

        # Split dataset into 80/20 training/test
        polygon_train = polygons.sample(frac=0.8, random_state=0)
        polygon_test = polygons.drop(polygon_train.index)

        train_features = polygon_train.copy()
        train_labels = train_features.pop('internal_nodes')

        test_features = polygon_test.copy()
        test_labels = test_features.pop('internal_nodes')

        polygon_model = tf.keras.Sequential([
            tf.keras.Input(shape=(13,)),
            layers.BatchNormalization(),
            layers.Dense(24, activation="relu"),

            layers.BatchNormalization(),
            layers.Dense(24, activation="relu"),

            layers.BatchNormalization(),
            layers.Dense(24, activation="relu"),

            layers.Dense(1)
        ])

        polygon_model.summary()
        polygon_model.compile(loss = tf.losses.MeanAbsoluteError(),
                              optimizer=tf.optimizers.Adam(
                                  # learning_rate=learning_rate, decay=weight_decay
                              ),
                              metrics=['accuracy']
                              )

        history = polygon_model.fit(train_features,
                                    train_labels,
                                    epochs=epochs,
                                    batch_size=batch_size,
                                    validation_data=(test_features, test_labels),
                                    verbose=2,
                                    )

        # Evaluate the model
        _, train_acc = polygon_model.evaluate(train_features, train_labels, verbose=0)
        _, test_acc = polygon_model.evaluate(test_features, test_labels, verbose=0)
        print('Train: %.3f, Test: %.3f' % (train_acc, test_acc))

        # plot history
        plt.plot(history.history['loss'], label='train')
        plt.plot(history.history['val_loss'], label='test')
        plt.legend()
        plt.show()

        # test_results = {}
        # test_results['polygon_model'] = polygon_model.evaluate(
        #     test_features,
        #     test_labels
        # )
        # print(test_results)

if __name__ == '__main__':
    start_time = time.time()
    nn1()
    print(f"Time elapsed: {round(time.time() - start_time, 5)} seconds")