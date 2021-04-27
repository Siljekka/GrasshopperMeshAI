import numpy as np


class GridPoint(x: float, y: float, score: float):
    def __init__(self):
        self.x = x
        self.y = y
        self.score = score

    def get_coordinates(self):
        return (self.x, self.y)

    def get_score(self):
        return self.score

    def set_score(self):
        self.score =


def point_grid(grid_resolution: int) -> list:
    # 1. Generate a point grid G of resolution N_g = n_g x n_g within a bounding box of size [-1.2, 1.2].
    #    Grid size 20 x 20 deemed adequate in papag. et al.
    # 2. Each grid point is assigned a score equalling the distance to the closest vertex.
    point_grid = []
    min_val = 1.2
    val_range = 2.4
    for i in range(grid_resolution+1):
        x_coord = i/grid_resolution * val_range - min_val
        for j in range(grid_resolution+1):
            y_coord = j/grid_resolution * val_range - min_val
            point_grid.append((x_coord, y_coord))

    point_grid
    return point_grid
