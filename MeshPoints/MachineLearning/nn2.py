import numpy as np
import pre_processing as pp
import matplotlib.path as mpl_path
from math import sqrt


class GridPoint():
    def __init__(self, pid: int, x: float, y: float, score: float = 0):
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
        # print(f"set grid point ({self._x}, {self._y}) neigbours to: {neighbours}")
        self._neighbours = neighbours

    @property
    def pid(self):
        return self._pid

    @property
    def score(self):
        return self._score

    @score.setter
    def score(self, score):
        print(f"set grid point ({self._x}, {self._y}) to score: {score}")
        self._score = score

    def __repr__(self):
        # return (f"{self._pid}")
        return (f"({round(self._x, 3)}, {round(self._y, 3)}, {self._score})")


def generate_point_grid(grid_resolution=20) -> list:
    VAL_RANGE = 2.64
    MIN_VAL = 1.32
    # 1. Generates a point grid G of resolution N_g = n_g x n_g within a bounding box of size [-1.2, 1.2].
    #    Grid size 20 x 20 deemed adequate in papag. et al. We add padding to help neighbour generation
    padded_grid = grid_resolution + 2
    point_grid = []
    point_id = 0
    for i in range(padded_grid, -1, -1):
        point_grid_row = []
        y_coord = i/padded_grid * VAL_RANGE - MIN_VAL
        for j in range(padded_grid+1):
            x_coord = j/padded_grid * VAL_RANGE - MIN_VAL
            point_grid_row.append(GridPoint(point_id, x_coord, y_coord))
            point_id += 1
        point_grid.append(point_grid_row)

    # Sets neighbours for each grid point (not including padding)
    for i in range(1, len(point_grid)-1):
        for j in range(1, len(point_grid)-1):
            neighbour_list = []
            neighbour_list.append(point_grid[i-1][j-1])     # bottom left
            neighbour_list.append(point_grid[i][j-1])       # bottom
            neighbour_list.append(point_grid[i+1][j-1])     # bottom right
            neighbour_list.append(point_grid[i+1][j])       # right
            neighbour_list.append(point_grid[i+1][j+1])     # top right
            neighbour_list.append(point_grid[i][j+1])       # top
            neighbour_list.append(point_grid[i-1][j+1])     # top left
            neighbour_list.append(point_grid[i-1][j])       # left
            point_grid[i][j].neighbours = neighbour_list

    # Flattens point grid list
    point_grid = [item for sublist in point_grid.copy() for item in sublist]

    # Removes padding from output by removing points that were not given neighbours
    point_grid = [x for x in point_grid.copy() if x.neighbours]

    return point_grid


def calculate_score(
        point_grid,
        internal_nodes,
        contour_points: np.array) -> None:
    # 0. Define contour (for finding points inside)
    contour_path = mpl_path.Path(contour_points)

    # For each point in point grid, iterate over internal nodes and calulcate euclidean distance.
    # If point is not in contour, set distance to 1000 (inf)
    for point in point_grid:
        score = 1000
        point_coordinates = point.get_coordinates()
        if contour_path.contains_point(point_coordinates):
            for internal_node in internal_nodes:
                distance = sqrt((point_coordinates[0]-internal_node[0])
                                ** 2 + (point_coordinates[1] - internal_node[1])**2)
                if distance < score:
                    score = distance
        point.score = score


if __name__ == "__main__":
    # np.array(pg).shape
    # print(pg)
    pg = generate_point_grid()
    contour = pp.create_random_ngon(6)
    random_points = pp.create_random_ngon(6)
    calculate_score(pg, random_points, contour)
    print(pg)
