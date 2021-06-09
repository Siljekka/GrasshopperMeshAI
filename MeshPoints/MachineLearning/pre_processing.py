import math
from math import cos, sin, pi
from random import random
import numpy as np
import matplotlib.pyplot as plt
import gmsh
import csv
import sys
import string


def simple_procrustes(contour: np.array) -> dict:
    n = contour.shape[0]  # number of nodes in input
    reference = create_regular_ngon(n)

    # Translate input surface, P*, to the origin.
    translation = np.mean(contour, 0)
    p_centered = contour - translation

    # Calculate the "centered Euclidean norm" of P*
    p_norm = np.linalg.norm(p_centered)

    # Scale centered_ P* by n/norm (Kabsch algorithm)
    scale_factor = math.sqrt(n) / p_norm
    p_scaled = p_centered * scale_factor

    # Apply "Singular Value Decomposition (SVD)" to A = P_centered.T . reference => UCV^T
    a = p_scaled.T @ reference
    u, c, vt = np.linalg.svd(a)

    # Get optimal rotation matrix, R = U*V^T
    rotation_matrix = u @ vt

    # Transformed contour given by P = p_scaled * rotation_matrix
    transformed_contour = p_scaled @ rotation_matrix

    assert np.allclose(
        contour,
        transformed_contour @ rotation_matrix.T * (1 / scale_factor) + translation,
    ), "Input contour does not match inverse transformed, transformed contour."

    return {
        "contour": transformed_contour,
        "scale": (1 / scale_factor),
        "rotation": rotation_matrix.T,
        "translation": translation,
    }


def create_regular_ngon(number_of_sides: int) -> np.array:
    """
    Creates a regular n-gon with circumradius = 1.
    This is equivalent to inscribing a regular polygon in a unit circle.
    """
    polygon = []
    for n in range(number_of_sides):
        polygon.append(
            [cos(2 * pi * n / number_of_sides), sin(2 * pi * n / number_of_sides)]
        )

    return np.array(polygon)


def create_random_ngon(number_of_sides: int) -> np.array:
    """
    Creates a random n-gon with a point placed in each sector of a unit disc divided into N sectors (number of sides).
    As per Papagiannopoulos et al. we do not want values of r too close to the origin.
    """

    polygon = []
    quantile = 1 / number_of_sides
    exclude = 0.2
    for n in range(number_of_sides):
        r = random() * (1 - exclude) + exclude  # mapping random from 0->1 to ex->1
        theta = 2 * pi * (n / number_of_sides + random() * quantile)
        polygon.append([r * cos(theta), r * sin(theta)])

    return np.array(polygon)


def create_random_displaced_ngon(number_of_sides: int) -> np.array:
    polygon = []
    quantile = 2 * pi / number_of_sides
    exclude = 0.2
    for n in range(number_of_sides):
        r = random() * (1 - exclude) + exclude  # mapping random from 0->1 to ex->1
        theta = 2 * pi * n / number_of_sides + random() * quantile
        polygon.append([r * cos(theta), r * sin(theta)])

    x_disp = random() * 1000 - 500
    y_disp = random() * 1000 - 500

    polygon = np.array(polygon) * random() * 100 + np.array([x_disp, y_disp])

    return polygon


def gmsh_settings() -> None:
    # Display settings
    gmsh.option.set_number("Mesh.Nodes", 1)
    gmsh.option.set_number("Mesh.NodeSize", 10)
    gmsh.option.set_number("Mesh.NodeLabels", 1)

    # Console settings
    gmsh.option.set_number("General.Verbosity", 0)


def mesh_contour(contour: np.array, target_edge_length: float):
    # Meshes a contour with a given edge length
    # Create a new model
    gmsh.model.add("1")

    for i, p in enumerate(contour):
        gmsh.model.geo.add_point(p[0], p[1], 0, target_edge_length, i)

    # Create curves
    curve_ids = []
    for j in range(len(contour)):
        if j == len(contour) - 1:
            gmsh.model.geo.add_line(j, 0, j)
        else:
            gmsh.model.geo.add_line(j, j + 1, j)
        # Constrain each edge curve to only have two points (endpoints)
        gmsh.model.geo.mesh.set_transfinite_curve(j, 2)
        curve_ids.append(j)

    # Add all ids of curves to define a curve loop
    gmsh.model.geo.add_curve_loop(curve_ids, 1)

    # Add curve loop with id=1 as a plane surface
    gmsh.model.geo.add_plane_surface([1], 1)

    # Synchronize CAD entities (point, line, surface) with the gmsh-model
    gmsh.model.geo.synchronize()
    gmsh.option.set_number("Mesh.CharacteristicLengthFactor", target_edge_length)
    gmsh.option.set_number("Mesh.Algorithm", 2)

    # Generate 2D mesh
    mesh = gmsh.model.mesh.generate(2)

    mesh_nodes = gmsh.model.mesh.get_nodes()
    num_nodes_total = mesh_nodes[0][-1]
    internal_nodes_count = int(num_nodes_total - contour.shape[0])

    # Get the coordinates of inserted internal nodes
    features = np.append(contour.copy(), internal_nodes_count)
    if internal_nodes_count > 0:
        internal_nodes_with_z = mesh_nodes[1][-internal_nodes_count * 3 :]
        internal_nodes = [x for i, x in enumerate(internal_nodes_with_z) if i % 3 != 2]
        features = np.append(features, internal_nodes)

    # GUI (buggy on macOS)
    # if "-nopopup" not in sys.argv:
    #     gmsh.fltk.run()
    gmsh.clear()

    return features


def plot_polygon(np_coords: np.array, style="") -> None:
    """
    See https://matplotlib.org/2.1.1/api/_as_gen/matplotlib.pyplot.plot.html for styles.
    """
    # Reshaping input and appending the first element to get a circular plot of the polygons
    coords = [[x[0] for x in np_coords], [y[1] for y in np_coords]]
    coords[0].append(coords[0][0])
    coords[1].append(coords[1][0])

    plt.plot(coords[0], coords[1], style, zorder=10)
    # Draw the first point as a red x
    # plt.plot(coords[0][0], coords[1][0], 'rx')


def generate_dataset(
    dataset_size: int, edge_count: int, target_edge_length: float
) -> list:
    # This helper function generates a dataset containing where each row has:
    #   - the coordinates of a contour
    #   - the number of nodes inserted by a reference mesher
    #   - the coordinates of said nodes

    # Start a gmsh API session
    gmsh.initialize()

    # suppress console output during generation
    gmsh_settings()
    dataset = []
    for i in range(dataset_size):
        print("meshing contour", i + 1, "of", dataset_size, end="\r")
        random_ngon = create_random_ngon(edge_count)
        transformed_polygon = simple_procrustes(random_ngon)

        dataset.append(mesh_contour(transformed_polygon["contour"], target_edge_length))

    # Close API call
    gmsh.finalize()
    print("\n")
    return dataset


def mesh_to_csv(features, dataset_size: int, number_of_sides: int) -> None:
    # Columns needed:
    #  contour nodes ( xi | yi) | target_edge_length | num inner nodes
    with open(
        f"data/{number_of_sides}-gon-mesh-with-internal-nodes-big.csv", "w", newline=""
    ) as file:
        writer = csv.writer(file)
        header = []
        for i in range(1, number_of_sides + 1):
            header.append(f"x{i}")
            header.append(f"y{i}")
        header.append("target_edge_length")
        header.append("internal_node_count")
        [header.append(x) for x in list(string.ascii_lowercase)]
        writer.writerow(header)

        print("\n")
        for i, contour in enumerate(features):
            print("writing", i, "of", dataset_size, end="\r")
            writer.writerow(contour)
