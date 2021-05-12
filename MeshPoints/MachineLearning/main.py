import pre_processing as pp
import dataset_generation as dg
import numpy as np
import matplotlib.pyplot as plt
import time
import gmsh

if __name__ == "__main__":
    # # ======== For plotting purposes ========
    # sample_polygon = pp.create_random_ngon(6)
    # p = pp.procrustes(sample_polygon)
    # pp.plot_polygon(p['transformed_contour'])
    # pp.plot_polygon(pp.create_regular_ngon(6))
    # plt.gca().set_aspect("equal", adjustable="box")
    # plt.show()

    # dataset_test = dg.generate_dataset(
    #     dataset_size, num_sides, target_edge_length)
    # dg.single_edge_length_mesh_to_csv(dataset_test, dataset_size, num_sides)

    # *** GMSH-tests ***
    # gmsh.initialize()
    # pp.gmsh_settings()
    # while True:
    #     sample_polygon = pp.create_random_ngon(num_sides)
    #     # pp.mesh_contour_quad(sample_polygon, target_edge_length)
    #     pp.mesh_contour_quad_all_lc(sample_polygon)
    # gmsh.finalize()

    # === Pre-processing ===
    num_sides = 6
    target_edge_length = 0.2
    # 1700 * 7 =~ 12_000 (wanted dataset size for 6-gons)
    dataset_size = 12_000
    # Trimeshing
    dataset_test = pp.generate_dataset_all_lc(dataset_size, num_sides)
    pp.mesh_to_csv(dataset_test, dataset_size, num_sides)
    # === Quad meshing ===
    # dataset = pp.generate_dataset_quad_all_lc(dataset_size, num_sides)
    # pp.quad_mesh_to_csv(dataset, dataset_size, num_sides)

    # === Neural networks ===
    # import nn_tests as nn
    # # import nn2
    # start_time = time.time()
    # nn.nn1()
    # # nn.nn1_single_edge_length(epochs=3000)
    # print(f"Time elapsed: {round(time.time() - start_time, 5)} seconds")
