import matplotlib.path as mpl_path
import numpy as np
import csv
from math import sqrt

"""
|---------------------------------------------------
|       CONSTANTS
|---------------------------------------------------
"""

GRID_RESOLUTION = 20
PATCH_SIZE = 2
VAL_RANGE = 2.86
GRID_SCORE_IF_OUTSIDE_CONTOUR = 1.5


"""
|---------------------------------------------------
|       CLASSES AND METHODS
|---------------------------------------------------
"""


class GridPoint:
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

    def mean_neighbourhood_score(self):
        # Average score of the grid point and all its neighbours
        neighbourhood_score = self._score
        if self._neighbours:
            for neighbour in self._neighbours:
                # Points on the edge have neighbours with score = 0 due to padding.
                # To exclude them from contention we give an arbitrary penalty.
                # Barring a statistical anomaly, no grid points should have score == 0.
                if neighbour._score == 0:
                    return 100
                else:
                    neighbourhood_score += neighbour._score
        return neighbourhood_score / (len(self._neighbours) + 1)

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
        # print(f"set grid point ({round(self._x, 3)}, {round(self._y, 3)}) to score: {round(score, 3)}")
        self._score = score

    def __repr__(self):
        # return (f"({round(self._x, 3)}, {round(self._y, 3)}, {self._score})")
        return f"{round(self._score, 3)}"


def generate_point_grid() -> list[list[GridPoint]]:
    """
    Output:
        - a point grid, G, of resolution 20x20 filled with GridPoint-objects with coordinates
          in the range ~[-1.4, -1.4] -> [1.4, 1.4].

    The grid is structured as a list of list where both the first row and column correspond
    to the GridPoint with the lowest value (and vice verse for the last row and column).

    Padding is added intermittently to simplify the generation of neighbours of a GridPoint.
    These padding points are removed before returning the list.
    """

    # Create list of coordinates from 0->VAL_RANGE, then shift it to be mirrored on 0.
    # Grid resolution + 2 for padding
    x_coordinates = (
        np.linspace(0.0, VAL_RANGE, num=GRID_RESOLUTION + 2, endpoint=True)
        - VAL_RANGE / 2
    )
    y_coordinates = (
        np.linspace(0.0, VAL_RANGE, num=GRID_RESOLUTION + 2, endpoint=True)
        - VAL_RANGE / 2
    )

    point_grid = []
    point_id = 0
    for y in y_coordinates:
        point_grid_row = []
        for x in x_coordinates:
            point_grid_row.append(GridPoint(point_id, x, y))
            point_id += 1
        point_grid.append(point_grid_row)

    # Sets neighbours for each grid point (not including padding)
    for y in range(1, len(point_grid) - 1):
        for x in range(1, len(point_grid) - 1):
            neighbour_list = []
            neighbour_list.append(point_grid[y - 1][x - 1])  # bottom left
            neighbour_list.append(point_grid[y - 1][x])  # bottom
            neighbour_list.append(point_grid[y - 1][x + 1])  # bottom right
            neighbour_list.append(point_grid[y][x + 1])  # right
            neighbour_list.append(point_grid[y + 1][x + 1])  # top right
            neighbour_list.append(point_grid[y + 1][x])  # top
            neighbour_list.append(point_grid[y + 1][x - 1])  # top left
            neighbour_list.append(point_grid[y][x - 1])  # left
            point_grid[y][x].neighbours = neighbour_list

    # Removes padding points from output by removing points that were not given neighbours
    cleaned_point_grid = []
    for row in point_grid:
        cleaned_row = [x for x in row if x.neighbours]
        # Completely empty rows (first and last) are also removed
        if cleaned_row:
            cleaned_point_grid.append(cleaned_row)

    return cleaned_point_grid


def calculate_score(
    point_grid: list[list[GridPoint]],
    internal_nodes: list[tuple[float]],
    contour_points: list[tuple[float]],
) -> None:
    """
    Input:
        - a point grid (list of lists)
        - coordinates of the internal nodes in the mesh (list of tuples)
        - coordinates of the contour vertices (list of tuples)

    For each GridPoint in the point grid, we iterate over the internal nodes and calculate
    the euclidean distances (scores) to each internal node. The shortest distance is the
    score of the GridPoint.

    If a GridPoint is not inside the contour, we set it's score to an arbitrary penalty.

    Current modifications (arbitrarly selected, but shown to give good prediction results):
        - points outside are mulitplied by 3.
        - points inside with score < 1 are raised to a power of 1.4
    """
    # Define contour (for finding points inside)
    contour_path = mpl_path.Path(contour_points)

    for row in point_grid:
        for point in row:
            score = 100
            point_coordinates = point.get_coordinates()

            for internal_node in internal_nodes:
                distance = sqrt(
                    (point_coordinates[0] - internal_node[0]) ** 2
                    + (point_coordinates[1] - internal_node[1]) ** 2
                )
                if distance < score:
                    score = distance
            if contour_path.contains_point(point_coordinates):
                if score < 1:
                    point.score = score ** 1.4
                else:
                    point.score = score
            else:
                point.score = 3 * score


def generate_patches(point_grid):
    """
    Generates patches of the given patch size (NxN) containing the GridPoint-objects
    of the patch.
    """

    patches = []
    for row in range(0, GRID_RESOLUTION, PATCH_SIZE):
        for col in range(0, GRID_RESOLUTION, PATCH_SIZE):
            p1 = point_grid[row][col]
            p2 = point_grid[row][col + 1]
            p3 = point_grid[row + 1][col]
            p4 = point_grid[row + 1][col + 1]

            patches.append([p1, p2, p3, p4])

    return patches


def generate_patch_training_dataset_from_contour(contour, point_grid):
    """This method processes the patches of a trained grid to turn the raw data
    of the point grid to patch data usable in training NN2. It returns a list of lists
    where each list contains all the necessary data for a row -- see output.

    Input:
        - coordinates of a contour
        - a point grid

    Output:
        - a list of list where each row of the list contains data for one patch:
            - the input contour coordinates [x1, y1 ... xn, yn],
            - the coordinates of the corners of the patch [x1, y1 ... x4, y4]
            - the grid scores of the corners of the patch [gs1 gs2 gs3 gs4]
    """

    patches = generate_patches(point_grid)

    dataset = []
    for patch in patches:
        patch_data = []

        patch_coordinates = [
            coordinate for point in patch for coordinate in point.get_coordinates()
        ]
        patch_grid_scores = [point.score for point in patch]

        patch_data.extend(contour)
        patch_data.extend(patch_coordinates)
        patch_data.extend(patch_grid_scores)

        dataset.append(patch_data)

    return dataset


def generate_patch_collection(
    dataset: pd.DataFrame, edge_count: int
) -> list[list[float]]:
    """
    Input:
        - a pandas.DataFrame where each row consists of:
            - the coordinates of a contour
            - the coordinates of any number of internal nodes

    Output:
        - a list of list where the rows contains all the "patch data" from all the given contours

    """

    patch_collection = []
    for raw_row in dataset.to_numpy():
        # Remove NaN values inserted during meshing
        row = raw_row[~np.isnan(raw_row)]

        # contour_coordinates [x1, y1 ... xn, yn] (n*2 values)
        # internal_node_coordinates [x1, y1 ... x_i, y_i] (i*2 values)
        contour_coordinates_flat_list = row[: edge_count * 2]
        internal_nodes_list = row[edge_count * 2 :]

        # Turn flat list into list of tuples (used in calculate_score)
        contour = list(
            zip(contour_coordinates_flat_list[::2], contour_coordinates_flat_list[1::2])
        )
        internal_nodes = list(zip(internal_nodes_list[::2], internal_nodes_list[1::2]))

        point_grid = generate_point_grid()
        calculate_score(point_grid, internal_nodes, contour)

        contour_patches = generate_patch_training_dataset_from_contour(
            contour_coordinates_flat_list, point_grid
        )

        patch_collection.extend(contour_patches)

    return patch_collection


def write_patch_collection_to_csv(
    patch_collection: list[list[float]], filename: str, edge_count: int
) -> None:
    with open(filename, "w", newline="") as file:
        writer = csv.writer(file)
        header = []
        for i in range(1, edge_count + 1):
            header.append(f"x{i}")
            header.append(f"y{i}")
        header.extend(
            [
                "gx1",
                "gy1",
                "gx2",
                "gy2",
                "gx3",
                "gy3",
                "gx4",
                "gy4",
                "sg1",
                "sg2",
                "sg3",
                "sg4",
            ]
        )
        writer.writerow(header)
        pc_len = len(patch_collection)
        for i, row in enumerate(patch_collection):
            print(f"writing patch {i+1} of {pc_len} to file {filename}", end="\r")
            writer.writerow(row)


def generate_patch_collection_and_write_to_csv(
    dataset, filename, edge_count=6, internal_count=2
) -> None:
    """Simple wrapper function for patch collection and csv-writing. Basically useless."""
    patch_collection = generate_patch_collection(dataset, edge_count)
    write_patch_collection_to_csv(patch_collection, filename, edge_count)


def generate_internal_nodes_from_grid_score(
    point_grid, target_internal_node_count: int
) -> list:
    """
    Input:
        - point grid with calculated/precited grid scores
        - target number of internal nodes within a contour,

    Output:
        - a list of the (x, y)-coordinates of the internal nodes interpolated from the grid scores.

    The interpolation procedure is as follows:
        1. Select the grid point with the lowest _neighbour_ score, meaning the sum of itself and
           all it's neighbours are the lowest in the grid.
        2. Give each point in the neighbourhood a weight based on its score: lower value -> higher weight.
           This is essentially a version of Inverse distance weighting (https://w.wiki/3TVw).
        3. Return the interpolated x- and y-coordinates of the internal node.
        4. Remove the entire neighbourhood + the neighbourhood's neighbours from being selected
           later.
    """
    internal_nodes = []
    # flattened point grid
    point_grid = [point for row in point_grid for point in row]

    while len(internal_nodes) < target_internal_node_count:

        min_grid_point = min(point_grid, key=lambda x: x.mean_neighbourhood_score())
        neighbours = min_grid_point.neighbours

        neighbourhood = [neighbour for neighbour in min_grid_point.neighbours]
        neighbourhood.append(min_grid_point)

        # Calculate weights of each node in the quadrant.
        total_score = sum(point.score ** 16 for point in neighbourhood)
        weights = np.array(
            [(total_score - point.score ** 16) / total_score for point in neighbourhood]
        )
        total_weight = sum(weights)

        # Interpolation
        interpolated_x = 0
        interpolated_y = 0
        for i, p in enumerate(neighbourhood):
            interpolated_x += weights[i] * p.x
            interpolated_y += weights[i] * p.y
        interpolated_x /= total_weight
        interpolated_y /= total_weight

        internal_nodes.append((interpolated_x, interpolated_y))

        # We remove the current minimum GridPoint and all its neighbours,
        # + all neighbours' neighbours from the point grid.
        if len(internal_nodes) != target_internal_node_count:
            points_to_exclude_by_id = {min_grid_point.pid}
            for neighbour in neighbours:
                neighbour_neighbours = neighbour.neighbours
                points_to_exclude_by_id.update(
                    [point.pid for point in neighbour_neighbours]
                )
                points_to_exclude_by_id.add(neighbour.pid)

            point_grid[:] = [
                point
                for point in point_grid
                if point.pid not in points_to_exclude_by_id
            ]

    return internal_nodes
