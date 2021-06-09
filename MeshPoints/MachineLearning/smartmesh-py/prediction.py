import tensorflow as tf
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mpl_path
import point_grid as pg
import pre_processing as pp
from math import sqrt

# ------------------------------------------
#            PREDICTION METHODS
# ------------------------------------------


def grid_prediction(contour: np.array, grid_model: tf.keras.Model) -> list[GridPoint]:
    """
    This method takes a contour and a grid model suited for predicitons on contours
    of the given size. An empty point grid is instantiated and filled with
    predictions, patch by patch.

    Returns a fully predicted point grid.
    """

    # Build point grid and patches
    point_grid = pg.generate_point_grid()
    patches = pg.generate_patches(point_grid)

    # Flatten contour list
    contour_list = [coordinate for point in contour for coordinate in point]

    for patch in patches:
        # Get patch coordinates
        patch_data = [
            coordinate for point in patch for coordinate in point.get_coordinates()
        ]

        # Prepare features for prediction and cast to correct data structure
        raw_features = np.append(contour_list, patch_data)
        features = np.expand_dims(raw_features, axis=0)

        # Prediction casted to numpy list
        prediction = grid_model(features).numpy()[0]

        # Write prediction to point grid
        for i, p in enumerate(patch):
            p.score = prediction[i]

    return point_grid


def direct_prediction(contour: np.array, direct_model: tf.keras.Model) -> list[float]:

    # Flatten contour list of tuples
    contour_list = [coordinate for point in contour for coordinate in point]

    features = np.expand_dims(contour_list, axis=0)
    prediction = direct_model(features).numpy()[0]

    return prediction


# ------------------------------------------
#            EVALUATION METHODS
# ------------------------------------------


def calc_distance_to_closest(
    contour: np.array, predicted_points: list, reference_points: list
) -> list[float]:
    """
    Calculates the minimal euclidean distance from each predicted
    internal point to the reference internal points.

    Assigns a negative value if point is outside of the contour.

    Returns a list of distances
    """
    error_list = []

    # Define contour (for finding points inside)
    contour_path = mpl_path.Path(contour)

    for p_point in predicted_points:
        # Set an arbitrary high distance
        error = 100

        # Loop over reference points to find the shortest distance
        for r_point in reference_points:
            shortest_distance = sqrt(
                (p_point[0] - r_point[0]) ** 2 + (p_point[1] - r_point[1]) ** 2
            )
            if shortest_distance < error:
                error = shortest_distance

        # Assign points outside contour a negative value for counting invalid points
        error = error if contour_path.contains_point(p_point) else -error
        error_list.append(error)

    return error_list


def calc_average_mean_and_worst_error(error_list: list[list[float]]) -> list[float]:
    """
    This method calculates the average over a sample of the mean and worst errors* of
    internal node predictions of meshes. Also counts the ratio of predicted nodes placed
    outside of it's respective contour.
    *this method is agnostic to the type of error input.

    Returns a list containing (Avg Mean Error, Avg Worst Error, Invalid point ratio)
    """
    sum_e_worst = 0  # sum of the largest errors of each mesh
    sum_e_mean = 0  # sum of the mse of each mesh
    outside = 0
    for mesh_error in error_list:
        for i, point_error in enumerate(mesh_error):
            # Negative values means predicted point was outside contour.
            # Count it and make it positive.
            if point_error < 0:
                outside += 1
                mesh_error[i] = -point_error

        sum_e_worst += max(mesh_error)
        sum_e_mean += sum(mesh_error) / len(mesh_error)

    e_worst = sum_e_worst / len(error_list)
    e_mean = sum_e_mean / len(error_list)

    return [round(e_mean, 3), round(e_worst, 3), outside / len(error_list)]


# Normalize contour with Procrustes and flatten list of tuples
# procrustes = pp.simple_procrustes(contour)
# transformed_contour = procrustes["contour"]
