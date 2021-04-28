import numpy as np


class GridPoint():

    def __init__(self, pid: int, x: float, y: float, score: float = 0):
        self._pid = pid
        self._x = x
        self._y = y
        self._score = score
        self._neighbours = []

    def coordinates(self):
        return (self._x, self._y)

    def set_neighbours(self, neighbours: list):
        self._neighbours = neighbours

    @property
    def pid(self):
        return self._pid

    @property
    def score(self):
        return self._score

    @score.setter
    def score(self, new_score):
        print(f"set grid point ({self._x}, {self._y}) to score: {new_score}")
        self._score = new_score

    def __repr__(self):
        # return (f"{self._pid}")
        return (f"({round(self._x, 3)}, {round(self._y, 3)}, {self._score})")


def generate_point_grid(grid_resolution=20) -> list:
    # 1. Generate a point grid G of resolution N_g = n_g x n_g within a bounding box of size [-1.2, 1.2].
    #    Grid size 20 x 20 deemed adequate in papag. et al. We add padding to help neighbour generation
    # 2. Generate neighbours of each point
    padded_grid = grid_resolution + 2
    point_grid = []
    min_val = 1.32
    val_range = 2.64
    point_id = 0
    for i in range(padded_grid+1, -1, -1):
        point_grid_row = []
        y_coord = i/padded_grid * val_range - min_val
        for j in range(padded_grid+1):
            x_coord = j/padded_grid * val_range - min_val
            point_grid_row.append(GridPoint(point_id, x_coord, y_coord))
            point_id += 1
        point_grid.append(point_grid_row)

    return point_grid


def set_neighbours_(point_grid):
    for i in range(1, len(point_grid)-1):

        for j in range(1, len(point_grid)-1):
            neighbour_list = []
            neighbour_list.append(point_grid[i-1, j-1].pid)  # bottom left
            neighbour_list.append(point_grid[i, j-1].pid)   # bottom
            neighbour_list.append(point_grid[i+1, j-1].pid)  # bottom right
            neighbour_list.append(point_grid[i+1, j].pid)   # right
            neighbour_list.append(point_grid[i+1, j+1].pid)  # top right
            neighbour_list.append(point_grid[i, j+1].pid)   # top
            neighbour_list.append(point_grid[i-1, j+1].pid)  # top left
            neighbour_list.append(point_grid[i-1, j].pid)  # left
            point_grid[i][j].set_neighbours(neighbour_list)


if __name__ == "__main__":
    pg = generate_point_grid()
    set_neighbours(pg)
