import pre_processing as pp
import dataset_generation as dg
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    # # ======== For plotting purposes ========
    # plt.gca().set_aspect("equal", adjustable="box")
    # plt.show()

    num_sides = 6
    target_edge_length = 0.4
    dataset_size = 1_333  # * 9 =~ 12_000
    # dataset_test = dg.generate_dataset(dataset_size, num_sides, target_edge_length)
    # dg.single_edge_length_mesh_to_csv(dataset_test, dataset_size, num_sides)

    dataset_test = dg.generate_dataset_all_lc(dataset_size, num_sides)
    dg.mesh_to_csv(dataset_test, dataset_size, num_sides)
