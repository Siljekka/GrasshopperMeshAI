import tensorflow as tf
from tensorflow.keras import layers
import pandas as pd
import numpy as np
import pre_processing as pp
import point_grid as pg

TARGET_EDGE_LENGTH = 0.4
DATASET_SIZE = {
    4: 10,
    6: 12000,
    8: 100,
    10: 48000,
    12: 95000,
}


def pre_processing(edge_count: int) -> pd.DataFrame:
    # Creates a pandas DataFrame for a given edge count
    dataset_size = DATASET_SIZE[edge_count]

    meshed_contours = pp.generate_dataset(
        dataset_size, edge_count, TARGET_EDGE_LENGTH
    )

    meshed_contours_dataframe = pd.DataFrame(meshed_contours)

    return meshed_contours_dataframe


def NN1_model_setup(edge_count: int) -> tf.keras.models.Model:
    model = tf.keras.Sequential([
        tf.keras.Input(shape=(edge_count*2,)),
        layers.BatchNormalization(),

        layers.Dense(edge_count*4, activation='relu'),
        layers.BatchNormalization(),

        layers.Dense(edge_count*4, activation='relu'),
        layers.BatchNormalization(),

        layers.Dense(1),
    ])

    return model


def NN1_training(edge_count: int, raw_data: pd.DataFrame) -> tf.keras.callbacks.History:
    # ===========================
    #          DATASET
    # ===========================

    dataset = raw_data.iloc[:, :edge_count*2+1]

    training_set = dataset.sample(frac=0.85)
    test_set = dataset.drop(training_set.index)

    train_features = training_set.iloc[:, :12]
    train_labels = training_set.iloc[:, -1:]

    test_features = test_set.iloc[:, :12]
    test_labels = test_set.iloc[:, -1:]

    # ===========================
    #         MODEL SETUP
    # ===========================
    epochs = 3000
    batch_size = 512

    model_path = f"auto-model/nn1-{edge_count}"
    lr_schedule = tf.keras.optimizers.schedules.ExponentialDecay(
        1e-2, 10_000, 3e-4
    )

    model = NN1_model_setup(edge_count)
    model.compile(
        loss=tf.losses.MeanAbsoluteError(),
        optimizer=tf.optimizers.Adam(
            learning_rate=lr_schedule,
        ),
    )

    # ===========================
    #          TRAINING
    # ===========================
    checkpoint = tf.keras.callbacks.ModelCheckpoint(
        model_path, monitor='val_loss', verbose=0, save_best_only=True, mode='min'
    )
    history = model.fit(train_features,
                        train_labels,
                        epochs=epochs,
                        batch_size=batch_size,
                        validation_split=0.18,
                        verbose=2,
                        callbacks=[checkpoint],
                        )

    # ===========================
    #        EVALUATION
    # ===========================
    train_acc = model.evaluate(
        train_features, train_labels, verbose=0)
    test_acc = model.evaluate(
        test_features, test_labels, verbose=0)

    print(f"=== Loss of nn1 with edge count {edge_count}")
    print(f"Training loss: {train_acc}, Test loss: {test_acc}")

    return history


def NN2_model_setup(edge_count: int) -> tf.keras.models.Model:
    model = tf.keras.Sequential([
        tf.keras.Input(shape=(edge_count*2+8,)),
        tf.keras.layers.BatchNormalization(),

        tf.keras.layers.Dense(edge_count*2+4, activation='relu'),
        tf.keras.layers.BatchNormalization(),

        tf.keras.layers.Dense(edge_count*2+4, activation='relu'),
        tf.keras.layers.BatchNormalization(),

        tf.keras.layers.Dense(4),
    ])

    return model


def NN2_training(edge_count: int, internal_node_count: int, raw_data: list):
    # ===========================
    #          DATASET
    # ===========================
    dataset = pd.DataFrame(raw_data)

    training_set = dataset.sample(frac=0.85)
    test_set = dataset.drop(training_set.index)

    # Split dataset into features and labels; last 4 (grid scores)
    train_features = training_set.iloc[:, :-4]
    train_labels = training_set.iloc[:, -4:]

    test_features = test_set.iloc[:, :-4]
    test_labels = test_set.iloc[:, -4:]
    print(test_features)


def NN2_wrapper(edge_count: int, raw_data: pd.DataFrame) -> list:
    model_histories = []

    # Loop through different internal node counts
    # for inc in range(1, 8):
    # INTERNAL NODE COUNT = INC
    inc = 1

    # ===========================
    #          DATASET
    # ===========================
    # Extract the internal node count column and use it to filter the dataset
    # to only contain the rows with inc == wanted inc
    inc_index = edge_count*2
    inc_df = raw_data.iloc[:, inc_index]
    df = raw_data[inc_df == float(inc)].dropna(axis=1, how="all")
    dataset = df.drop(inc_index, axis=1)

    patch_dataset = pg.generate_patch_collection(
        dataset,
        edge_count=edge_count,
        internal_count=inc
    )

    history = NN2_training(edge_count, inc, patch_dataset)
    # return model_histories
    pass


if __name__ == "__main__":
    edge_count = 4

    # 1. Create dataset
    print(f"=== Creating dataset for edge count = {edge_count} ===\n")
    mesh_data = pre_processing(edge_count)
    # 2. Train internal node count on nn1_training
    nn1_training_history = []
    nn2_training_history = []
    print(NN2_wrapper(edge_count, mesh_data))

    # 3. Create patch dataset for internal node count. Limit to 1-7 for 8gon
