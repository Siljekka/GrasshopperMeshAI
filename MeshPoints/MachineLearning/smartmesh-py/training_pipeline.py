import tensorflow as tf
from tensorflow.keras import layers
from datetime import datetime
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pre_processing as pp
import point_grid as pg
import os

# Suppress tensorflow runtime messages in terminal
# Doesn't work.
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
tf.get_logger().setLevel("ERROR")

TARGET_EDGE_LENGTH = 0.4
DATASET_SIZE = {
    4: 10,
    6: 12000,
    8: 24000,
    10: 48000,
    12: 95000,
}


def create_meshed_contour_dataset(edge_count: int) -> pd.DataFrame:
    """Creates the base meshed contour dataset used for training NN1 and NN2."""
    dataset_size = DATASET_SIZE[edge_count]

    meshed_contours = pp.generate_dataset(dataset_size, edge_count, TARGET_EDGE_LENGTH)

    meshed_contours_dataframe = pd.DataFrame(meshed_contours)

    return meshed_contours_dataframe


def NN1_model_setup(edge_count: int) -> tf.keras.models.Model:
    model = tf.keras.Sequential(
        [
            tf.keras.Input(shape=(edge_count * 2,)),
            layers.BatchNormalization(),
            layers.Dense(edge_count * 4, activation="relu"),
            layers.BatchNormalization(),
            layers.Dense(edge_count * 4, activation="relu"),
            layers.BatchNormalization(),
            layers.Dense(1),
        ]
    )

    return model


def NN1_training(edge_count: int, raw_data: pd.DataFrame) -> tf.keras.callbacks.History:
    # ===========================
    #          DATASET
    # ===========================

    dataset = raw_data.iloc[:, : edge_count * 2 + 1]

    training_set = dataset.sample(frac=0.85)
    test_set = dataset.drop(training_set.index)

    train_features = training_set.iloc[:, : edge_count * 2]
    train_labels = training_set.iloc[:, -1:]

    test_features = test_set.iloc[:, : edge_count * 2]
    test_labels = test_set.iloc[:, -1:]

    # ===========================
    #         MODEL SETUP
    # ===========================
    epochs = 3000
    batch_size = 512

    model_path = f"model/nn1-{edge_count}gon"

    model = NN1_model_setup(edge_count)

    lr_schedule = tf.keras.optimizers.schedules.PolynomialDecay(
        1e-2, epochs * batch_size // 2, 3e-4
    )
    model.compile(
        loss=tf.losses.MeanAbsoluteError(),
        optimizer=tf.optimizers.Adam(
            learning_rate=lr_schedule,
        ),
    )

    # ===========================
    #          TRAINING
    # ===========================

    log_dir = f"logs/{model_path}-" + datetime.now().strftime("%Y%m%d-%H%M%S")
    tensorboard = tf.keras.callbacks.TensorBoard(log_dir=log_dir, histogram_freq=1)

    checkpoint = tf.keras.callbacks.ModelCheckpoint(
        model_path, monitor="val_loss", verbose=0, save_best_only=True, mode="min"
    )

    history = model.fit(
        train_features,
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
    train_acc = model.evaluate(train_features, train_labels, verbose=0)
    test_acc = model.evaluate(test_features, test_labels, verbose=0)

    print(f"=== Loss of NN1 with edge count {edge_count} ===")
    print(f"Training loss: {train_acc}, Test loss: {test_acc} \n")

    return history


def NN2_model_setup(edge_count: int) -> tf.keras.models.Model:
    model = tf.keras.Sequential(
        [
            tf.keras.Input(shape=(edge_count * 2 + 8,)),
            tf.keras.layers.BatchNormalization(),
            tf.keras.layers.Dense(edge_count * 2 + 4, activation="relu"),
            tf.keras.layers.BatchNormalization(),
            tf.keras.layers.Dense(edge_count * 2 + 4, activation="relu"),
            tf.keras.layers.BatchNormalization(),
            tf.keras.layers.Dense(4),
        ]
    )

    return model


def NN2_training(edge_count: int, raw_data: list) -> tf.keras.callbacks.History:
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

    model_path = f"model/nn2-{edge_count}gon"

    model = NN2_model_setup(edge_count)

    lr_schedule = tf.keras.optimizers.schedules.PolynomialDecay(
        1e-1, epochs * batch_size // 2, 3e-4
    )
    model.compile(
        loss=tf.losses.MeanSquaredError(),
        optimizer=tf.optimizers.Adam(
            learning_rate=lr_schedule,
        ),
    )

    # ===========================
    #          TRAINING
    # ===========================

    checkpoint = tf.keras.callbacks.ModelCheckpoint(
        model_path, monitor="val_loss", verbose=0, save_best_only=True, mode="min"
    )
    early_stopping = tf.keras.callbacks.EarlyStopping(
        monitor="val_loss", mode="min", patience=epochs // 5, min_delta=0.001
    )

    log_dir = f"logs/{model_path}-" + datetime.now().strftime("%Y%m%d-%H%M%S")
    tensorboard = tf.keras.callbacks.TensorBoard(log_dir=log_dir, histogram_freq=1)

    history = model.fit(
        train_features,
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
    train_acc = model.evaluate(train_features, train_labels, verbose=0)
    test_acc = model.evaluate(test_features, test_labels, verbose=0)

    print(f"=== Loss of NN2 with edge count: {edge_count}===")
    print(f"Training loss: {train_acc}, Test loss: {test_acc}\n")

    return history


def NN2_wrapper(edge_count: int, raw_data: pd.DataFrame) -> tf.keras.callbacks.History:
    """
    This function creates the patch dataset for NN2 calls the NN2_training function.

    Returns a History object for the trained model.
    """

    # Index of the internal node count in raw_data
    inc_index = edge_count * 2

    # ===========================
    #    PATCH DATA GENERATION
    # ===========================
    print(f"=== Generating patch dataset for NN2 with edge count: {edge_count} ===\n")
    inc_df = raw_data.iloc[:, inc_index]

    # Drop rows with 0 internal nodes
    df = raw_data[inc_df != 0]

    # Drop column containing internal node count
    dataset = df.drop(inc_index, axis=1)

    patch_dataset = pg.generate_patch_collection(
        dataset,
        edge_count=edge_count,
    )

    print(f"=== Training NN2 for edge count: {edge_count} ===\n")
    history = NN2_training(edge_count, patch_dataset)

    return history


if __name__ == "__main__":
    edge_count = 6

    # 1. Create dataset
    print(f"=== Creating dataset for edge count: {edge_count} ===\n")
    mesh_data = create_meshed_contour_dataset(edge_count)

    # 2. Train NN1
    nn1_training_history = []
    print(f"=== Training NN1 for edge count: {edge_count} ===\n")
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
