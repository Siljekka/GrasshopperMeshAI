import pre_processing as pp
import dataset_generation as dg
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

    # pp.plot_polygon(procrustes["transformed_contour"], "ko-")
    # pp.plot_polygon(pp.create_regular_ngon(num_sides), "ro-")
    # pp.plot_polygon(test, ":y")

    # # ======== For plotting purposes ========
    # plt.gca().set_aspect("equal", adjustable="box")
    # plt.show()
    #pp.to_csv(procrustes, num_sides)

    num_sides = 6
    target_edge_length = 0.4
    dataset_size = 24_000
    dataset_test = dg.generate_dataset(dataset_size, num_sides, target_edge_length)
    dg.mesh_to_csv(dataset_test, dataset_size, num_sides)
