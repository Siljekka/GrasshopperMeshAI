import pre_processing as pp
import dataset_generation as dg
import numpy as np
import matplotlib.pyplot as plt
import time
import nn_tests as nn

if __name__ == "__main__":
    # # ======== For plotting purposes ========
    # plt.gca().set_aspect("equal", adjustable="box")
    # plt.show()

    # === Pre-processing ===
    num_sides = 6
    target_edge_length = 0.4
    dataset_size = 1_700  # 1700 * 7 =~ 12_000 (wanted dataset size for 6-gons)
    # # dataset_test = dg.generate_dataset(dataset_size, num_sides, target_edge_length)
    # # dg.single_edge_length_mesh_to_csv(dataset_test, dataset_size, num_sides)

    # Try removing bottom and max vals of target edge length
    # dataset_test = pp.generate_dataset_all_lc(dataset_size, num_sides)
    # pp.mesh_to_csv(dataset_test, dataset_size, num_sides)

    # print(pp.point_grid(20))

    # === Neural networks ===
    start_time = time.time()
    nn.nn1()
    print(f"Time elapsed: {round(time.time() - start_time, 5)} seconds")
