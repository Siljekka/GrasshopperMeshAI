import pre_processing as ct
import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    # test_contour = np.array(
    #     [[7.06, 0.55],
    #      [9.16, -0.99],
    #      [9.62, -3.25],
    #      [14.2, -0.95],
    #      [14.36, 3.31],
    #      [10.46, 3.53]])

    num_sides = 10
    test = ct.create_random_ngon(num_sides)
    procrustes = ct.procrustes(test)
    ct.plot_polygon(procrustes['transformed_contour'], 'ko-')
    ct.plot_polygon(ct.create_regular_ngon(num_sides), 'ko-')
    ct.plot_polygon(test, ':y')

    # ======== For plotting purposes ========
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()
