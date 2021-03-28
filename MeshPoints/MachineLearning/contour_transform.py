from math import cos, sin, pi
import numpy as np
from scipy import linalg, spatial
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
    p_mean = np.array([np.sum(contour[0]), np.sum(contour[1])])/n
    # i think this is ~0 always, so we could remove it
    q_mean = [np.sum(reg_ngon[0])/n, np.sum(reg_ngon[1])/n]

    # Calculate the "centered Euclidean norms" of P and Q
    # todo: is this a scalar or a vector (of scalars) ???
    p_norm = []
    q_norm = []
    for i in range(n):
        p_norm.append((contour[0][i] - p_mean[0])**2 +
                      (contour[1][i] - p_mean[1])**2)
        q_norm.append((reg_ngon[0][i])**2 + (reg_ngon[1][i])**2)

    p_norm = sum(p_norm)
    q_norm = sum(q_norm)

    # Scale P* and Q by the norms
    # p_scaled=np.zeros((2,n))
    # q_scaled=np.zeros((2,n))
    # for i in range(n):
    #     p_scaled[0][i] = contour[0][i]/p_norm
    #     p_scaled[1][i] = contour[1][i]/p_norm
    #     q_scaled[0][i] = reg_ngon[0][i]/q_norm
    #     q_scaled[1][i] = reg_ngon[1][i]/q_norm

    p_scaled = contour/p_norm
    q_scaled = reg_ngon/q_norm
    # q_scaled = np.array([reg_ngon[i]/p_norm[i] for i in range(n)])
    # q_scaled = np.array([reg_ngon[0]/q_norm, reg_ngon[1]/q_norm])

    # Apply "Singular Value Decomposition (SVD)" to A = q_scaled.T * p_scaled => A = UCV^T
    A = reg_ngon.copy().T @ contour
    U, C, Vt = linalg.svd(A)  # svd[0] = U, svd[1] = C, svd[2] = V^T

    # Get optimal rotation matrix, R, and scaling factor, S, of input contour.
    # S = q_norm * trace(C) where trace(C) = sum(C) because of linalg.svd outputting only non-zero values of C
    scaling_factor = q_norm / np.sum(C)

    # R = U*V^T
    rotation_matrix = U @ Vt

    # The coordinates of the transformed contour are given by P_trans = scaling_factor * p_scaled * rotation + q_mean
    # We add the q_mean values separately
    tmp = (scaling_factor * p_scaled) @ rotation_matrix
    tc_x = tmp[0] + q_mean[0]
    tc_y = tmp[1] + q_mean[1]
    transformed_contour = np.array([tc_x, tc_y])
    plot_ngon(q_scaled)
    plot_ngon(p_scaled)
    return {"transformed_contour": transformed_contour, "scaling_factor": scaling_factor}


def procrustes_superimposition_scipy(contour: np.array):
    n = contour.shape[1]
    reg_ngon = create_regular_ngon(n)

    # have to transform list of points to [x1, y1] ... [xn, yn]-shape
    # not [x1... xn], [y1 ... yn]

    reshaped_c =[]
    reshaped_ngon = []
    for i in range(n):
        reshaped_c.append([contour[0][i], contour[1][i]])
        reshaped_ngon.append([reg_ngon[0][i], reg_ngon[1][i]])


    transformed_contour = spatial.procrustes(reshaped_ngon, reshaped_c)

    return transformed_contour


def plot_ngon(np_coords: np.array):
    coords = np_coords.tolist()
    coords[0].append(coords[0][0])
    coords[1].append(coords[1][0])

    plt.plot(coords[0], coords[1])
    plt.gca().set_aspect('equal', adjustable='box')
