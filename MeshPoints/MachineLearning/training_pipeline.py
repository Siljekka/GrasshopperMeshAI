import tensorflow as tf
from tensorflow.keras import layers
from datetime import datetime
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pre_processing as pp
import point_grid as pg

TARGET_EDGE_LENGTH = 0.4
DATASET_SIZE = {
    4: 10,
    6: 12000,
    8: 24000,
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

    train_features = training_set.iloc[:, :edge_count*2]
    train_labels = training_set.iloc[:, -1:]

    test_features = test_set.iloc[:, :edge_count*2]
    test_labels = test_set.iloc[:, -1:]

    # ===========================
    #         MODEL SETUP
    # ===========================
    epochs = 3000
    batch_size = 512

    model_path = f"auto-model/nn1-{edge_count}gon"

    model = NN1_model_setup(edge_count)

    lr_schedule = tf.keras.optimizers.schedules.ExponentialDecay(
        1e-2, 10_000, 3e-4)
    model.compile(
        loss=tf.losses.MeanAbsoluteError(),
        optimizer=tf.optimizers.Adam(
            learning_rate=lr_schedule,
        ),
    )

    # ===========================
    #          TRAINING
    # ===========================
    print(f"=== Training NN1 for edge count: {edge_count} ===\n")

    log_dir = f"logs/{model_path}-" + datetime.now().strftime("%Y%m%d-%H%M%S")
    tensorboard = tf.keras.callbacks.TensorBoard(
        log_dir=log_dir, histogram_freq=1)

    checkpoint = tf.keras.callbacks.ModelCheckpoint(
        model_path, monitor='val_loss', verbose=0, save_best_only=True, mode='min')

    history = model.fit(train_features,
                        train_labels,
                        epochs=epochs,
                        batch_size=batch_size,
                        validation_split=0.18,
                        verbose=0,
                        callbacks=[checkpoint, tensorboard],
                        )

    # ===========================
    #        EVALUATION
    # ===========================
    train_acc = model.evaluate(
        train_features, train_labels, verbose=0)
    test_acc = model.evaluate(
        test_features, test_labels, verbose=0)

    print(f"=== Loss of NN1 with edge count {edge_count} ===")
    print(f"Training loss: {train_acc}, Test loss: {test_acc} \n")

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

    train_features = training_set.iloc[:, :-4]
    train_labels = training_set.iloc[:, -4:]

    test_features = test_set.iloc[:, :-4]
    test_labels = test_set.iloc[:, -4:]

    # ===========================
    #         MODEL SETUP
    # ===========================
    epochs = 5000
    batch_size = 512

    model_path = f"auto-model/nn2-{edge_count}gon-{internal_node_count}int"

    model = NN2_model_setup(edge_count)

    lr_schedule = tf.keras.optimizers.schedules.ExponentialDecay(
        1e-1, 10000, 1e-3)
    model.compile(loss=tf.losses.MeanSquaredError(),
                  optimizer=tf.optimizers.Adam(
                  learning_rate=lr_schedule,
                  ),
                  )

    # ===========================
    #          TRAINING
    # ===========================

    checkpoint = tf.keras.callbacks.ModelCheckpoint(
        model_path, monitor='val_loss', verbose=0, save_best_only=True, mode='min')
    early_stopping = tf.keras.callbacks.EarlyStopping(
        monitor='val_loss', mode='min', patience=epochs//5, min_delta=0.001)

    log_dir = f"logs/{model_path}-" + datetime.now().strftime("%Y%m%d-%H%M%S")
    tensorboard = tf.keras.callbacks.TensorBoard(
        log_dir=log_dir, histogram_freq=1)

    history = model.fit(train_features,
                        train_labels,
                        epochs=epochs,
                        batch_size=batch_size,
                        validation_split=0.18,
                        verbose=0,
                        callbacks=[checkpoint, early_stopping, tensorboard],
                        )

    # ===========================
    #        EVALUATION
    # ===========================
    train_acc = model.evaluate(
        train_features, train_labels, verbose=0)
    test_acc = model.evaluate(
        test_features, test_labels, verbose=0)

    print(
        f"=== Loss of NN2 with edge count {edge_count}, internal node count {internal_node_count} ===")
    print(f"Training loss: {train_acc}, Test loss: {test_acc}\n")

    return history


def NN2_wrapper(edge_count: int, raw_data: pd.DataFrame) -> list:
    """
    This function both creates the patch dataset for NN2 and trains NN2 on
    several different number of internal nodes.

    Returns a History object for each trained model.
    """

    model_histories = []
    inc_index = edge_count*2
    max_inc = raw_data.iloc[:, inc_index].max()
    print(max_inc)
    # Loop through different internal node counts
    for inc in range(1, int(max_inc)):
        print(
            f"=== Training NN2 for edge count: {edge_count} and inc: {inc} ===\n")
        # INTERNAL NODE COUNT = INC

        # ===========================
        #    PATCH DATA GENERATION
        # ===========================
        try:
            inc_df = raw_data.iloc[:, inc_index]
            df = raw_data[inc_df == float(inc)].dropna(axis=1, how="all")
            dataset = df.drop(inc_index, axis=1)

            patch_dataset = pg.generate_patch_collection(
                dataset,
                edge_count=edge_count,
                internal_count=inc
            )

            model_histories.append(NN2_training(
                edge_count, inc, patch_dataset))
        except KeyError:
            print(f"No examples of contours with {inc} internal nodes.")
            break
    # return model_histories
    return model_histories


if __name__ == "__main__":
    # Suppress tensorflow runtime messages in terminal
    import os
    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
    tf.get_logger().setLevel('ERROR')

    edge_count = 10

    # 1. Create dataset
    print(f"=== Creating dataset for edge count: {edge_count} ===\n")
    mesh_data = pre_processing(edge_count)

    # 2. Train NN1
    nn1_training_history = []
    nn1_training_history.append(NN1_training(edge_count, mesh_data))

    # 3. TRAIN NN2
    nn2_training_history = []
    nn2_training_history.append(NN2_wrapper(edge_count, mesh_data))

    # ======================
    #        PLOTTING
    #    (w/o tensorboard)
    # ======================
    # for history in nn1_training_history:
    #     plt.plot(history.history['loss'], label='training loss')
    #     plt.plot(history.history['val_loss'], label='validation loss')
    #     plt.legend()
    #     plt.xlabel("Epochs")
    #     plt.ylabel("MAE")

    # for history_set in nn2_training_history:
    #     for history in history_set:
    #         plt.plot(history.history['loss'], label='training loss')
    #         plt.plot(history.history['val_loss'], label='validation loss')
    #         plt.legend()
    #         plt.xlabel("Epochs")
    #         plt.ylabel("MSE")
