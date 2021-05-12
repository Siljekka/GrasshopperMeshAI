from math import cos, sin, pi
from random import random
import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
import gmsh
import csv
import sys
import string


def procrustes(contour: np.array) -> dict:
    """
    In the paper "How to teach neural networks to mesh: Application to 2-D simplical contours" [1]
    the authors specify a pre-processing step using Procrustes superimposition [2] to assist pattern
    recognition in the neural networks.

    The input contour, P*, with n sides is transformed to match a reference regular ngon, Q.
    The scaling factor says how much to scale the input contour to get the transformed contour

    [1] Papagiannopoulos, A.; Clausen, P., Avellan, F. (2020).
        "How to teach neural networks to mesh: Application to 2-D simplical contours"
    [2] Gower, J. C. (1975). "Generalized procrustes analysis".
    """
    n = contour.shape[0]  # number of nodes (and therefore sides) in input
    reference = create_regular_ngon(n)

    # Translate P* to the origin.
    # As Q is defined in the origin, we do not have to translate the reference polygon
    p_centered = contour - np.mean(contour, 0)

    # Calculate the "centered Euclidean norms" of P* and Q
    p_norm = np.linalg.norm(p_centered)
    q_norm = np.linalg.norm(reference)

    # Scale P* and Q by the norms
    p_scaled = p_centered / p_norm
    q_scaled = reference / q_norm

    # Apply "Singular Value Decomposition (SVD)" to A = q_scaled.T * p_scaled => A = UCV^T
    a = q_scaled.T @ p_scaled
    u, c, vt = linalg.svd(a)

    # Get the sigma, S = q_norm * trace(C); trace(C) = sum(C) in scipy.linalg.svd. Aka. scaling factor.
    sigma = q_norm * np.sum(c)

    # Get optimal rotation matrix, R = U*V^T
    rotation_matrix = u @ vt

    # The coordinates of the transformed contour are given by P_trans = scaling_factor * p_scaled * rotation
    transformed_contour = sigma * p_scaled @ rotation_matrix

    # The scale difference between the transformed and the original contour is s / p_norm;
    assert np.allclose(contour, p_norm * p_scaled + np.mean(contour, 0))
    scale = sigma / p_norm

    return {"transformed_contour": transformed_contour, "scale": scale}


def create_regular_ngon(number_of_sides: int) -> np.array:
    """
    Creates a regular n-gon with circumradius = 1.
    This is equivalent to inscribing a regular polygon in a unit circle.
    """
    polygon = []
    for n in range(number_of_sides):
        polygon.append(
            [cos(2 * pi * n / number_of_sides),
             sin(2 * pi * n / number_of_sides)]
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
        theta = 2 * pi * n / number_of_sides + random() * quantile
        polygon.append([r * cos(theta), r * sin(theta)])

    return np.array(polygon)


def gmsh_settings() -> None:
    # Display settings
    gmsh.option.set_number("Mesh.Nodes", 1)
    gmsh.option.set_number("Mesh.NodeSize", 10)
    gmsh.option.set_number("Mesh.NodeLabels", 1)

    # Console settings
    gmsh.option.set_number("General.Verbosity", 0)


def mesh_contour(contour: np.array, target_edge_length: float):
    # Meshes a contour with a given edge length
    gmsh.initialize()
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
    # CharacteristicLengthFactor
    gmsh.option.set_number(
        "Mesh.CharacteristicLengthFactor", target_edge_length)
    # gmsh.model.geo.mesh.setSize()
    gmsh.option.set_number("Mesh.Algorithm", 5)

    # Generate 2D mesh
    mesh = gmsh.model.mesh.generate(2)

    mesh_nodes = gmsh.model.mesh.get_nodes()
    num_nodes_total = mesh_nodes[0][-1]
    internal_nodes_count = int(num_nodes_total - contour.shape[0])

    features = np.append(contour.copy(), target_edge_length)
    features = np.append(features, internal_nodes_count)
    if internal_nodes_count > 0:
        internal_nodes_with_z = mesh_nodes[1][-internal_nodes_count*3:]
        internal_nodes = [x for i, x in enumerate(
            internal_nodes_with_z) if i % 3 != 2]
        features = np.append(features, internal_nodes)

    # GUI (buggy on macOS)
    # Display options:
    gmsh.option.set_number("Mesh.Nodes", 1)
    gmsh.option.set_number("Mesh.NodeSize", 10)
    gmsh.option.set_number("Mesh.NodeLabels", 1)
    # if "-nopopup" not in sys.argv:
    #     gmsh.fltk.run()

    gmsh.clear()
    gmsh.finalize()

    return features


def mesh_contour_quad(contour: np.array, target_edge_length: float):
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
    cl = gmsh.model.geo.add_curve_loop(curve_ids)

    # Add curve loop with id=1 as a plane surface
    pl = gmsh.model.geo.add_plane_surface([cl])

    # Synchronize CAD entities (point, line, surface) with the gmsh-model
    gmsh.model.geo.synchronize()
    gmsh.option.set_number(
        "Mesh.CharacteristicLengthFactor", target_edge_length)
    gmsh.option.set_number("Mesh.Algorithm", 5)

    # To generate quadrangles instead of triangles
    gmsh.model.mesh.setRecombine(2, pl)

    # Generate 2D mesh
    mesh = gmsh.model.mesh.generate(2)

    mesh_nodes = gmsh.model.mesh.get_nodes()
    num_nodes_total = mesh_nodes[0][-1]
    internal_nodes_count = int(num_nodes_total - contour.shape[0])

    features = np.append(contour.copy(), target_edge_length)
    features = np.append(features, internal_nodes_count)
    if internal_nodes_count > 0:
        internal_nodes_with_z = mesh_nodes[1][-internal_nodes_count*3:]
        internal_nodes = [x for i, x in enumerate(
            internal_nodes_with_z) if i % 3 != 2]
        features = np.append(features, internal_nodes)

    # GUI (buggy on macOS)
    # Display options:
    if "-nopopup" not in sys.argv:
        gmsh.fltk.run()

    gmsh.clear()

    return features


def mesh_contour_all_lc(contour: np.array):
    # Meshes one contour with all the edge lengths

    edge_lengths = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    features = []
    for k, target_edge_length in enumerate(edge_lengths):
        tmp_features = np.zeros(len(contour)+1)
        # Create a new model
        gmsh.model.add(f"{k}")

        for i, p in enumerate(contour):
            # (x, y, z, edge_length, id)
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

        # ! Unsure about this line !!!!
        gmsh.option.set_number(
            "Mesh.CharacteristicLengthFactor", target_edge_length)
        gmsh.option.set_number("Mesh.Algorithm", 2)

        # Generate 2D mesh
        mesh = gmsh.model.mesh.generate(2)

        # Extract features
        mesh_nodes = gmsh.model.mesh.get_nodes()
        num_nodes_total = mesh_nodes[0][-1]
        internal_nodes_count = int(num_nodes_total - contour.shape[0])

        # For NN1: contour nodes, target edge length, and num of internal nodes
        tmp_features = np.append(contour.copy(), target_edge_length)
        tmp_features = np.append(tmp_features, internal_nodes_count)
        if internal_nodes_count > 0:
            internal_nodes_with_z = mesh_nodes[1][-internal_nodes_count*3:]
            internal_nodes = [x for i, x in enumerate(
                internal_nodes_with_z) if i % 3 != 2]
            tmp_features = np.append(tmp_features, internal_nodes)

        # GUI (buggy on macOS)
        # Display options:
        # gmsh.option.set_number("Mesh.Nodes", 1)
        # gmsh.option.set_number("Mesh.NodeSize", 10)
        # gmsh.option.set_number("Mesh.NodeLabels", 1)
        # if "-nopopup" not in sys.argv:
        #     gmsh.fltk.run()

        # print(f"edge length: {target_edge_length}")
        # print(f"internal nodes: {internal_nodes}")

        # Write to file
        # gmsh.write("data/1.msh")
        features.append(tmp_features)

        # Clear model data
        gmsh.clear()

    return features


def mesh_contour_quad_all_lc(contour: np.array):
    # Meshes one contour with all the edge lengths
    edge_lengths = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    features = []
    for k, target_edge_length in enumerate(edge_lengths):
        tmp_features = np.zeros(len(contour)+1)
        # Create a new model
        gmsh.model.add(f"{k}")

        for i, p in enumerate(contour):
            # (x, y, z, edge_length, id)
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
        cl = gmsh.model.geo.add_curve_loop(curve_ids, 1)

        # Add curve loop with id=1 as a plane surface
        pl = gmsh.model.geo.add_plane_surface([cl])

        # Mesh settings
        gmsh.option.set_number(
            "Mesh.CharacteristicLengthFactor", target_edge_length)
        # gmsh.option.set_number("Mesh.Algorithm", 8)

        # Synchronize CAD entities (point, line, surface) with the gmsh-model
        gmsh.model.geo.synchronize()

        # To generate quadrangles instead of triangles
        gmsh.model.mesh.setRecombine(2, pl)

        # Generate 2D mesh
        mesh = gmsh.model.mesh.generate(2)

        # Extract features
        mesh_nodes = gmsh.model.mesh.get_nodes()
        num_nodes_total = mesh_nodes[0][-1]
        internal_nodes_count = int(num_nodes_total - contour.shape[0])
        # For NN1: contour nodes, target edge length, and num of internal nodes
        tmp_features = np.append(contour.copy(), target_edge_length)
        tmp_features = np.append(tmp_features, internal_nodes_count)
        if internal_nodes_count > 0:
            internal_nodes_with_z = mesh_nodes[1][-internal_nodes_count*3:]
            internal_nodes = [x for i, x in enumerate(
                internal_nodes_with_z) if i % 3 != 2]
            tmp_features = np.append(tmp_features, internal_nodes)

        # GUI (buggy on macOS)
        # Display options:
        # if "-nopopup" not in sys.argv:
        #     gmsh.fltk.run()

        # Write to file
        # gmsh.write("data/1.msh")
        features.append(tmp_features)

        # Clear model data
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

    plt.plot(coords[0], coords[1], style)


def generate_dataset(dataset_size: int, num_sides: int, target_edge_length: float) -> list:

    # Start a gmsh API session
    gmsh.initialize()

    # suppress console output during generation
    gmsh.option.set_number("General.Verbosity", 0)
    dataset = []
    for _ in range(dataset_size):
        test_polygon = create_random_ngon(num_sides)
        transformed_polygon = procrustes(test_polygon)

        dataset.append(mesh_contour(
            transformed_polygon["transformed_contour"], target_edge_length))

    # Close API call
    gmsh.finalize()

    return dataset


def generate_dataset_all_lc(dataset_size: int, num_sides: int) -> list:

    # Start a gmsh API session
    gmsh.initialize()
    gmsh_settings()
    dataset = []
    for i in range(dataset_size):
        print("meshing contour", i, "of", dataset_size, end="\r")
        test_polygon = create_random_ngon(num_sides)
        transformed_polygon = procrustes(test_polygon)

        dataset_nested = mesh_contour_all_lc(
            transformed_polygon["transformed_contour"])
        for l in dataset_nested:
            dataset.append(l)

    # Close API call
    gmsh.finalize()

    return dataset


def generate_dataset_quad_all_lc(dataset_size: int, num_sides: int) -> list:

    # Start a gmsh API session
    gmsh.initialize()
    gmsh_settings()
    dataset = []
    for i in range(dataset_size):
        print("meshing contour", i, "of", dataset_size, end="\r")
        test_polygon = create_random_ngon(num_sides)
        transformed_polygon = procrustes(test_polygon)

        dataset_nested = mesh_contour_quad_all_lc(
            transformed_polygon["transformed_contour"])
        for l in dataset_nested:
            dataset.append(l)

    # Close API call
    gmsh.finalize()

    return dataset


def contour_to_csv(triangles: np.array, number_of_sides: int) -> None:
    with open(f"data/{number_of_sides}-gon-dataset", "w", newline="") as file:
        writer = csv.writer(file)
        header = ["id"]
        for i in range(1, number_of_sides + 1):
            header.append(f"x{i}")
            header.append(f"y{i}")
        header.append("scaling_ratio")

        writer.writerow(header)


def single_edge_length_mesh_to_csv(features, dataset_size: int, number_of_sides: int) -> None:
    # Columns needed:
    #  contour nodes ( xi | yi) | target_edge_length | num inner nodes
    with open(f"data/{number_of_sides}-gon-lc-04-mesh-dataset.csv", "w", newline="") as file:
        writer = csv.writer(file)
        header = []
        for i in range(1, number_of_sides + 1):
            header.append(f"x{i}")
            header.append(f"y{i}")
        header.append("target_edge_length")
        header.append("internal_nodes")

        writer.writerow(header)

        for i, contour in enumerate(features):
            print("writing", i, "of", dataset_size, end="\r")
            writer.writerow(contour)


def mesh_to_csv(features, dataset_size: int, number_of_sides: int) -> None:
    # Columns needed:
    #  contour nodes ( xi | yi) | target_edge_length | num inner nodes
    with open(f"data/{number_of_sides}-gon-mesh-with-internal-nodes.csv", "w", newline="") as file:
        writer = csv.writer(file)
        header = []
        for i in range(1, number_of_sides + 1):
            header.append(f"x{i}")
            header.append(f"y{i}")
        header.append("target_edge_length")
        header.append("internal_node_count")
        [header.append(x) for x in list(string.ascii_lowercase)]
        writer.writerow(header)

        for i, contour in enumerate(features):
            print("writing", i, "of", dataset_size, end="\r")
            writer.writerow(contour)


def quad_mesh_to_csv(features, dataset_size: int, number_of_sides: int) -> None:
    # Columns needed:
    #  contour nodes ( xi | yi) | target_edge_length | num inner nodes
    with open(f"data/{number_of_sides}-gon-quad-mesh-dataset.csv", "w", newline="") as file:
        writer = csv.writer(file)
        header = []
        for i in range(1, number_of_sides + 1):
            header.append(f"x{i}")
            header.append(f"y{i}")
        header.append("target_edge_length")
        header.append("internal_node_count")
        [header.append(x) for x in list(string.ascii_lowercase)]
        writer.writerow(header)

        for i, contour in enumerate(features):
            print("writing", i, "of", dataset_size, end="\r")
            writer.writerow(contour)
