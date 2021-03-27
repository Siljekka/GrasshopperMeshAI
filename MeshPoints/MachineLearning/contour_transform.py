from math import cos, sin, pi
import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt


def create_regular_ngon(number_of_sides: int):
    """
    Creates a regular n-gon with circumradius = 1.
    This is equivalent to inscribing a regular polygon in a unit circle.
    """
    x = []
    y = []
    for n in range(number_of_sides):
        x.append(cos(2*pi*n/number_of_sides))
        y.append(sin(2*pi*n/number_of_sides))

    return np.array([x, y])


def create_random_ngon(number_of_sides: int):


    pass


def procrustes_superimposition(contour: np.array):
    n = contour.shape[1]  # number of nodes (and therefore sides) in input
    reg_ngon = create_regular_ngon(n)

    # Calculate the mean coordinate values of input contour, P, and regular polygon, Q
    p_mean = [np.divide(np.sum(contour[0]), n), np.divide(np.sum(contour[1]), n)]
    q_mean = [np.divide(np.sum(reg_ngon[0]), n), np.divide(np.sum(reg_ngon[1]), n)]

    # Calculate the "centered Euclidean norms" of P and Q
    p_norm = np.zeros((2,1))
    q_norm = np.zeros((2,1))
    for i in range(2):
        p_tmp = p_mean[i]
        q_tmp = q_mean[i]
        for j in range(n):
            p_norm[i] += (contour[i][j] - p_tmp)**2
            q_norm[i] += (reg_ngon[i][j] - q_tmp)**2

    # Scale P* and Q on the norms
    p_scaled = np.array([np.divide(contour[0], p_norm[0]), np.divide(contour[1], p_norm[1])])
    q_scaled = np.array([np.divide(reg_ngon[0], q_norm[0]), np.divide(reg_ngon[1], q_norm[1])])

    # Apply "Singular Value Decomposition (SVD)" to A = q_scaled.transpose * p_scaled => A = UCV
    a = np.matmul(q_scaled.transpose(), p_scaled)
    svd_decomp = linalg.svd(a)  # svd[0] = U, svd[1] = C, svd[2] = V

    # Get optimal rotation matrix, R, and scaling factor, S, of input contour.
    # S = q_norm * trace(C) where trace(C) = sum(C) because of linalg.svd outputting only non-zero values of C
    scaling_factor = q_norm[0] * np.sum(svd_decomp[1])

    # R = U*V
    rotation_matrix = np.matmul(svd_decomp[0], svd_decomp[2])

    # The coordinates of the transformed contour are given by P_trans = scaling_factor * p_scaled * rotation + q_mean
    # We add the q_mean values separately
    tmp =  np.matmul(scaling_factor * p_scaled, rotation_matrix)
    tc_x = np.add(tmp[0], q_mean[0])
    tc_y = np.add(tmp[1], q_mean[1])
    transformed_contour = np.array([tc_x, tc_y])
    plot_ngon(q_scaled)
    return transformed_contour


def plot_ngon(np_coords: np.array):
    coords = np_coords.tolist()
    coords[0].append(coords[0][0])
    coords[1].append(coords[1][0])

    plt.plot(coords[0], coords[1])
    plt.gca().set_aspect('equal', adjustable='box')
