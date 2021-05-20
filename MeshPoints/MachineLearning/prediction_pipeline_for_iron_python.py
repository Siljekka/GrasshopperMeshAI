from math import sqrt
import numpy as np
import matplotlib.path as mpl_path
import matplotlib.pyplot as plt
import pre_processing as pp
import csv
import pandas as pd
import tensorflow as tf


GRID_RESOLUTION = 20
PATCH_SIZE = 2
VAL_RANGE = 2.86
GRID_SCORE_IF_OUTSIDE_CONTOUR = 1.5


"""
|---------------------------------------------------
|       CLASSES AND METHODS
|---------------------------------------------------
"""


class GridPoint():
    def __init__(self, pid, x, y, score=0):
        self._pid = pid
        self._x = x
        self._y = y
        self._score = score
        self._neighbours = []

    def get_coordinates(self):
        return (self._x, self._y)

    @property
    def neighbours(self):
        return self._neighbours

    @neighbours.setter
    def neighbours(self, neighbours: list):
        self._neighbours = neighbours

    def mean_neighbourhood_score(self):
        # Average score of the grid point and all its neighbours
        neighbourhood_score = self._score
        if self._neighbours:
            for neighbour in self._neighbours:
                # Points on the edge have neighbours with score = 0 due to padding
                if neighbour._score == 0:
                    return 100
                else:
                    neighbourhood_score += neighbour._score
        return neighbourhood_score/(len(self._neighbours) + 1)

    @property
    def x(self):
        return self._x

    @property
    def y(self):
        return self._y

    @property
    def pid(self):
        return self._pid

    @property
    def score(self):
        return self._score

    @score.setter
    def score(self, score):
        self._score = score


def generate_patches(point_grid):
    patches = []
    for row in range(0, GRID_RESOLUTION, PATCH_SIZE):
        for col in range(0, GRID_RESOLUTION, PATCH_SIZE):
            p1 = point_grid[row][col]
            p2 = point_grid[row][col+1]
            p3 = point_grid[row+1][col]
            p4 = point_grid[row+1][col+1]

            patches.append([p1, p2, p3, p4])

    return patches


def generate_point_grid():
    x_coordinates = np.linspace(
        0.0, VAL_RANGE, num=GRID_RESOLUTION+2, endpoint=True) - VAL_RANGE/2
    y_coordinates = np.linspace(
        0.0, VAL_RANGE, num=GRID_RESOLUTION+2, endpoint=True) - VAL_RANGE/2

    point_grid = []
    point_id = 0
    for y in y_coordinates:
        point_grid_row = []
        for x in x_coordinates:
            point_grid_row.append(GridPoint(point_id, x, y))
            point_id += 1
        point_grid.append(point_grid_row)

    # Sets neighbours for each grid point (not including padding)
    for y in range(1, len(point_grid)-1):
        for x in range(1, len(point_grid)-1):
            neighbour_list = []
            neighbour_list.append(point_grid[y-1][x-1])     # bottom left
            neighbour_list.append(point_grid[y-1][x])       # bottom
            neighbour_list.append(point_grid[y-1][x+1])     # bottom right
            neighbour_list.append(point_grid[y][x+1])       # right
            neighbour_list.append(point_grid[y+1][x+1])     # top right
            neighbour_list.append(point_grid[y+1][x])       # top
            neighbour_list.append(point_grid[y+1][x-1])     # top left
            neighbour_list.append(point_grid[y][x-1])       # left
            point_grid[y][x].neighbours = neighbour_list

    # Removes padding points from output by removing points that were not given neighbours
    cleaned_point_grid = []
    for row in point_grid:
        cleaned_row = [x for x in row if x.neighbours]
        # Completely empty rows (first and last) are also removed
        if cleaned_row:
            cleaned_point_grid.append(cleaned_row)

    return cleaned_point_grid


def generate_internal_nodes_from_grid_score(point_grid, target_internal_node_count: int) -> list:

    internal_nodes = []
    # flattened point grid
    point_grid = [point for row in point_grid for point in row]

    while len(internal_nodes) < target_internal_node_count:

        min_grid_point = min(point_grid, key=lambda x: x.score)
        neighbours = min_grid_point.neighbours

        gridpoint_quadrant = [
            neighbour for neighbour in min_grid_point.neighbours]
        gridpoint_quadrant.append(min_grid_point)

        # Calculate weights of each node in the quadrant.
        total_score = sum(point.score**16 for point in gridpoint_quadrant)
        weights = np.array([(total_score-point.score**16) /
                           total_score for point in gridpoint_quadrant])
        total_weight = sum(weights)

        # Interpolation
        interpolated_x = 0
        interpolated_y = 0
        for i, p in enumerate(gridpoint_quadrant):
            interpolated_x += weights[i]*p.x
            interpolated_y += weights[i]*p.y
        interpolated_x /= total_weight
        interpolated_y /= total_weight

        internal_nodes.append((interpolated_x, interpolated_y))

        # We remove the current minimum GridPoint and all its neighbours, and all neighbours' neighbours from the point grid.
        if len(internal_nodes) != target_internal_node_count:
            points_to_exclude_by_id = {min_grid_point.pid}
            for neighbour in neighbours:
                neighbour_neighbours = neighbour.neighbours
                points_to_exclude_by_id.update(
                    [point.pid for point in neighbour_neighbours])
                points_to_exclude_by_id.add(neighbour.pid)

            point_grid[:] = [
                point for point in point_grid if point.pid not in points_to_exclude_by_id]

    return internal_nodes


# Pre-processing method for nn2.
def prediction_pipeline_nn2(contour, internal_node_count):

    # Build point grid
    point_grid = generate_point_grid()
    patches = generate_patches(point_grid)

    # Pre-process contour with procrustes superimposition
    transformed_contour = pp.procrustes(contour)['transformed_contour']
    transformed_contour_data = [
        coordinate for point in transformed_contour for coordinate in point]

    # Load model from file
    model = tf.keras.models.load_model('model/grid-score-zerofour-lc-2')
    # model.summary()

    # Get patch data
    for patch in patches:
        patch_data = [
            coordinate for point in patch for coordinate in point.get_coordinates()]

        # Define prediction data
        features = np.append(transformed_contour_data, patch_data)
        prediction_data = np.expand_dims(features, axis=0)

        # Predict
        prediction_result = model(prediction_data).numpy()

        # Write predicted score to the point grid
        for i, p in enumerate(patch):
            p.score = prediction_result[0][i]

    return point_grid


if __name__ == "__main__":
    internal_node_count = 2
    contour = np.array([[1, 1], [4, 4], [0, 6], [-4, 3], [-4, 1], [0, 0]])

    predicted_point_grid = prediction_pipeline_nn2(
        contour, internal_node_count)
    pg_internal_nodes = generate_internal_nodes_from_grid_score(
        predicted_point_grid, internal_node_count)
    
    print(pg_internal_nodes)
    
