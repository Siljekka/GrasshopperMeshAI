from math import cos, sin, pi
from random import random
import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
import gmsh
import sys
import csv


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
        theta = 2 * pi * n / number_of_sides + random() * quantile
        polygon.append([r * cos(theta), r * sin(theta)])

    return np.array(polygon)


def mesh_contour(contour: np.array, target_edge_length: float):
    """
    Mesh the transformed contour using Gmsh [1] as the reference mesher.

    We do not place points along the edges of the contour. It is crucial
    that the input points are defined

    [1] C. Geuzaine and J.-F. Remacle. Gmsh: a three-dimensional finite element
        mesh generator with built-in pre- and post-processing facilities. 2009.

    output: .msh-file containing mesh data
    """
    # Start a gmsh API session
    gmsh.initialize()

    # Create a new model
    gmsh.model.add("1")

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

    # Generate 2D mesh
    gmsh.model.mesh.generate(2)

    # GUI (buggy on macOS)
    if "-nopopup" not in sys.argv:
        gmsh.fltk.run()

    # Write to file
    gmsh.write("data/1.msh")

    # Close API call
    gmsh.finalize()


def plot_polygon(np_coords: np.array, style="") -> None:
    """
    See https://matplotlib.org/2.1.1/api/_as_gen/matplotlib.pyplot.plot.html for styles.
    """
    # Reshaping input and appending the first element to get a circular plot of the polygons
    coords = [[x[0] for x in np_coords], [y[1] for y in np_coords]]
    coords[0].append(coords[0][0])
    coords[1].append(coords[1][0])

    plt.plot(coords[0], coords[1], style)
