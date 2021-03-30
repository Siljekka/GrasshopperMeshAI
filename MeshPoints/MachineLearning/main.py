import pre_processing as pp
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    # test_contour = np.array(
    #     [[7.06, 0.55],
    #      [9.16, -0.99],
    #      [9.62, -3.25],
    #      [14.2, -0.95],
    #      [14.36, 3.31],
    #      [10.46, 3.53]])

    num_sides = 5
    test = pp.create_random_ngon(num_sides)
    procrustes = pp.procrustes(test)
    # pp.plot_polygon(procrustes["transformed_contour"], "ko-")
    # pp.plot_polygon(pp.create_regular_ngon(num_sides), "ko-")
    # pp.plot_polygon(test, ":y")

    # # ======== For plotting purposes ========
    # plt.gca().set_aspepp("equal", adjustable="box")
    # plt.show()
    pp.to_csv(procrustes, num_sides)
    # pp.mesh_contour(procrustes["transformed_contour"], 1.0)
