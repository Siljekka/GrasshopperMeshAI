from math import cos, sin, pi
from random import random
import numpy as np
from scipy import linalg, spatial
import matplotlib.pyplot as plt


def procrustes(contour: np.array) -> dict:
    """
    In the paper "How to teach neural networks to mesh: Application to 2-D simplical contours" [1]
    the authors specify a pre-processing step using Procrustes superimposition [2] to assist pattern
    recognition in the neural networks.

    The input contour, P*, with n sides is transformed to match a reference regular ngon, Q.

    [1] Papagiannopoulos, A.; Clausen, P., Avellan, F. (2020).
        "How to teach neural networks to mesh: Application to 2-D simplical contours"
    [2] Gower, J. C. (1975). "Generalized procrustes analysis".
    """
    n = contour.shape[0]  # number of nodes (and therefore sides) in input
    reference = create_regular_ngon(n)

    # Translate P* to the origin.
    p_mean = contour - np.mean(contour, 0)
    assert(sum(np.mean(reference, 0)) < 0.000001)  # As Q is defined in the origin, we do not have to move it
    q_mean = reference - np.mean(reference, 0)

    # Calculate the "centered Euclidean norms" of P* and Q
    p_norm = np.linalg.norm(p_mean)
    q_norm = np.linalg.norm(q_mean)

    # Scale P* and Q by the norms
    p_scaled = p_mean/p_norm
    q_scaled = q_mean/q_norm

    # Apply "Singular Value Decomposition (SVD)" to A = q_scaled.T * p_scaled => A = UCV^T
    A = q_scaled.T @ p_scaled
    U, C, Vt = linalg.svd(A)

    # Get scaling factor, S = q_norm * trace(C); trace(C) = sum(C) in scipy.linalg.svd.
    scaling_factor = q_norm * np.sum(C)

    # Get optimal rotation matrix, R = U*V^T
    rotation_matrix = U @ Vt

    # The coordinates of the transformed contour are given by P_trans = scaling_factor * p_scaled * rotation + ~q
    transformed_contour = scaling_factor * p_scaled @ rotation_matrix

    return {"transformed_contour": transformed_contour, "scaling_factor": scaling_factor}


def create_regular_ngon(number_of_sides: int) -> np.array:
    """
    Creates a regular n-gon with circumradius = 1.
    This is equivalent to inscribing a regular polygon in a unit circle.
    """
    polygon = []
    for n in range(number_of_sides):
        polygon.append([
            cos(2 * pi * n / number_of_sides),
            sin(2 * pi * n / number_of_sides)
        ])

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
        polygon.append([
            r * cos(theta),
            r * sin(theta)
        ])

    return np.array(polygon)

def plot_polygon(np_coords: np.array) -> None:
    # Reshaping input and appending the first element to get a circular plot of the polygons
    coords = [[x[0] for x in np_coords], [y[1] for y in np_coords]]
    coords[0].append(coords[0][0])
    coords[1].append(coords[1][0])

    plt.plot(coords[0], coords[1])

