from math import sqrt
import numpy as np
import matplotlib.path as mpl_path
import csv


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
        #print(f"set grid point ({round(self._x, 3)}, {round(self._y, 3)}) to score: {round(score, 3)}")
        self._score = score

    def __repr__(self):
        # return (f"{self._pid}")
        return (f"({round(self._x, 3)}, {round(self._y, 3)}, {self._score})")


def generate_point_grid():
    """
    |   Output: 
    |       - a point grid, G, of resolution 20x20 filled with GridPoint-objects with coordinates
    |         in the range ~[-1.3, -1.3] -> [1.3, 1.3]. 
    |        
    |   The grid is structured as a list of list where both the first row and column correspond         
    |   to the GridPoint with the lowest value (and vice verse for the last row and column).
    |
    |   Padding is added intermittently to simplify the generation of neighbours of a GridPoint.     
    |   These padding points are removed before returning the list.

    !!! function is hardcoded to 20x20 in the range ~[-1.3, -1.3] -> [1.3, 1.3] to stop myself from messing up !!!
    """

    VAL_RANGE = 2.86
    GRID_RESOLUTION = 20

    # Create list of coordinates from 0->VAL_RANGE, then shift it to be mirrored on 0.
    # Grid resolution + 2 for padding
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
        if cleaned_row is not None:
            cleaned_point_grid.append(cleaned_row)

    return cleaned_point_grid


def calculate_score(point_grid, internal_nodes, contour_points) -> None:
    """
    |   Input: 
    |       - a point grid (list of lists) 
    |       - coordinates of the internal nodes in the mesh (list of tuples)
    |       - coordinates of the contour vertices (list of tuples)
    |
    |   For each GridPoint in the point grid, we iterate over the internal nodes and calculate
    |   the euclidean distances (scores) to each internal node. The shortest distance is the
    |   score of the GridPoint. 
    |
    |   If a GridPoint is not inside the contour, we set it's score = 2. 
    |
    """
    # Define contour (for finding points inside)
    contour_path = mpl_path.Path(contour_points)

    for row in point_grid:
        for point in row:
            score = 2
            point_coordinates = point.get_coordinates()
            if contour_path.contains_point(point_coordinates):
                for internal_node in internal_nodes:
                    distance = sqrt((point_coordinates[0]-internal_node[0])
                                    ** 2 + (point_coordinates[1] - internal_node[1])**2)
                    if distance < score:
                        score = distance
            point.score = score


def generate_patch_training_dataset_from_contour(contour, point_grid):
    """
    |   Input:
    |       - coordinates of a contour
    |       - a point grid
    |
    |   Output: 
    |       - a list of list where each row of the list contains data for one patch:
    |           - the input contour coordinates [x1 y1 ... xn yn], 
    |           - the coordinates of the corners of the patch [x1 y1 ... x4 y4]
    |           - the grid scores of the corners of the patch [gs1 gs2 gs3 gs4]
    """

    patch_size = 2
    grid_resolution = 20

    patches = []
    for row in range(0, grid_resolution, patch_size):
        for col in range(0, grid_resolution, patch_size):
            p1 = point_grid[row][col]
            p2 = point_grid[row][col+1]
            p3 = point_grid[row+1][col]
            p4 = point_grid[row+1][col+1]

            patch_coordinates = [p1.x, p1.y, p2.x, p2.y,
                                 p3.x, p3.y, p4.x, p4.y]

            patch = contour
            patch.extend(patch_coordinates)
            patch.extend([p1.score, p2.score, p3.score, p4.score])

            patches.append(patch)

    return patches


def generate_patch_collection(dataset):
    """
    !!! only functional for 6gons with 2 internal nodes !!!
    |   Input:
    |       - a pandas.DataFrame where each row consists of:
    |           - the coordinates of a contour
    |           - the coordinates of a set number of internal nodes
    |
    |   Output: 
    |       - a list of list where the rows contains all the "patch data" from all the given contours
    |
    |   Used for csv generation; feeds into write_patch_collection_to_csv.
    |
    """

    patch_collection = []
    for row in dataset.values.tolist():
        # For a 6-gon with 2 internal nodes, a row is structured like:
        # contour_coordinates [x1 -> y6] (12 values)
        # internal_node_count = 2 (1 value)
        # internal_node_coordinates [x1 -> y2] (4 values)
        contour_coordinates_flat_list = row[:12]
        internal_nodes_list = row[-4:]

        # Turn flat list into list of tuples (used in calculate_score)
        contour = list(
            zip(contour_coordinates_flat_list[::2], contour_coordinates_flat_list[1::2]))
        internal_nodes = list(
            zip(internal_nodes_list[::2], internal_nodes_list[1::2]))

        point_grid = generate_point_grid()
        calculate_score(point_grid, internal_nodes, contour)

        contour_patches = generate_patch_training_dataset_from_contour(
            contour_coordinates_flat_list, point_grid)

        patch_collection.extend(contour_patches)

    return patch_collection

# Method for writing a patch collection to csv


def write_patch_collection_to_csv(patch_collection):
    with open(f"data/patches-dataset.csv", "w", newline="") as file:
        writer = csv.writer(file)
        header = ['x1', 'y1',
                  'x2', 'y2',
                  'x3', 'y3',
                  'x4', 'y4',
                  'x5', 'y5',
                  'x6', 'y6',
                  'gx1', 'gy1',
                  'gx2', 'gy2',
                  'gx3', 'gy3',
                  'gx4', 'gy4',
                  'sg1', 'sg2',
                  'sg3', 'sg4'
                  ]
        writer.writerow(header)
        for row in patch_collection:
            writer.writerow(row)


def generate_internal_nodes_from_grid_score(point_grid, target_internal_node_count: int) -> list:
    """
    |   Input:
    |       - point grid with calculated/precited grid scores
    |       - target number of internal nodes within a contour,
    | 
    |   Output:
    |       - a list of the (x, y)-coordinates of the internal nodes interpolated from the grid scores.
    |
    |-------------------------------
    |   Interpolation
    |-------------------------------
    |   The interpolation process is based on two ideas. Given a grid point with minimal score, we assume;
    |    1. The actual point lies within the quadrant that contains the diagonal neighbour with the lowest score.
    |    2. The location of the actual point within the quadrant can be approximated using a procedure similar
    |       to how one would find the center of mass in mechanics. We weight each node in the quadrant using the 
    |       inverse square of their grid scores.
    """
    internal_nodes = []
    # flattened point grid
    point_grid = [point for row in point_grid for point in row]

    # The quadrants-dictionary contains the indices of all the points in a
    # quadrant given the index of a "diagonal" neighbour to the GridPoint-object.
    # 6  5  4
    # 7  P  3
    # 0  1  2
    quadrants = {
        0: [7, 0, 1],
        2: [1, 2, 3],
        4: [3, 4, 5],
        6: [5, 6, 7],
    }

    while len(internal_nodes) < target_internal_node_count:

        min_grid_point = min(point_grid, key=lambda x: x.score)
        neighbours = min_grid_point.neighbours

        # Decides which quadrant the point is located in.
        lowest_score_diagonal_node = neighbours[0]
        lowest_score_diagonal_node_index = 0
        for i in [2, 4, 6]:
            if neighbours[i].score < lowest_score_diagonal_node.score:
                lowest_score_diagonal_node = neighbours[i]
                lowest_score_diagonal_node_index = i

        # Generates a list containing the grid points in the quadrant.
        gridpoint_quadrant = [neighbours[i]
                              for i in quadrants[lowest_score_diagonal_node_index]]
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

        # We remove the current minimum GridPoint and all its neighbours from the point grid.
        if len(internal_nodes) != target_internal_node_count:
            points_to_exclude_by_id = {min_grid_point.pid}
            points_to_exclude_by_id.update([point.pid for point in neighbours])

            point_grid[:] = [
                point for point in point_grid if point.pid not in points_to_exclude_by_id]

    return internal_nodes
