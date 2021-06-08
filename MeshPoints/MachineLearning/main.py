import pre_processing as pp
import numpy as np
import matplotlib.pyplot as plt
import time
import gmsh
# import tensorflow as tf

if __name__ == "__main__":
    # # ======== For plotting purposes ========
    # while True:
    #     sample_polygon = pp.create_random_ngon(6)
    #     p = pp.procrustes(sample_polygon)
    #     pp.plot_polygon(p['transformed_contour'])
    #     pp.plot_polygon(pp.create_regular_ngon(6))
    #     plt.gca().set_aspect("equal", adjustable="box")
    #     plt.grid(True)
    #     plt.show()

    # p = pp.create_random_displaced_ngon(6)
    # pp.plot_polygon(p)
    # plt.show()

    # =========================================================
    # PROCRUSTES DEMO
    y1, y2, y3, y4, y5, y6 = [], [], [], [], [], []
    x1, x2, x3, x4, x5, x6 = [], [], [], [], [], []
    for _ in range(1000):
        pro = pp.create_random_ngon(6)
        pro = pp.simple_procrustes(pro)['contour']
        x1.append(pro[0][0])
        y1.append(pro[0][1])
        x2.append(pro[1][0])
        y2.append(pro[1][1])
        x3.append(pro[2][0])
        y3.append(pro[2][1])
        x4.append(pro[3][0])
        y4.append(pro[3][1])
        x5.append(pro[4][0])
        y5.append(pro[4][1])
        x6.append(pro[5][0])
        y6.append(pro[5][1])

    plt.scatter(x1, y1, c='r', marker='D', alpha=0.5)
    plt.scatter(x2, y2, c='b', marker='D', alpha=0.5)
    plt.scatter(x3, y3, c='g', marker='D', alpha=0.5)
    plt.scatter(x4, y4, c='r', marker='D', alpha=0.5)
    plt.scatter(x5, y5, c='b', marker='D', alpha=0.5)
    plt.scatter(x6, y6, c='g', marker='D', alpha=0.5)

    plt.grid(True)
    plt.gca().set_aspect("equal", adjustable="box")
    params = {'axes.axisbelow': True,
              'xtick.labelsize': 40,
              'ytick.labelsize': 40,

              }
    plt.rcParams.update(params)
    plt.show()
    # ========================================================

    # *** GMSH-tests ***
    # num_sides = 8
    # gmsh.initialize()
    # pp.gmsh_settings()
    # while True:
    #     sample_polygon = pp.create_random_ngon(num_sides)
    #     pro = pp.procrustes(sample_polygon)['transformed_contour']
    #     # pp.mesh_contour_quad(sample_polygon, target_edge_length)
    #     pp.mesh_contour(sample_polygon, 0.4)
    # gmsh.finalize()

    # === Pre-processing ===
    # num_sides = 8
    # target_edge_length = 0.4
    # # 1700 * 7 =~ 12_000 (wanted dataset size for 6-gons)
    # dataset_size = 24_000
    # dataset_test = pp.generate_dataset(
    #     dataset_size, num_sides, target_edge_length)
    # pp.single_edge_length_mesh_to_csv(dataset_test, dataset_size, num_sides)

    # Trimeshing
    # dataset_test = pp.generate_dataset_all_lc(dataset_size, num_sides)
    # pp.mesh_to_csv(dataset_test, dataset_size, num_sides)

    # === Neural networks ===
    # import nn_tests as nn
    # # import nn2
    # start_time = time.time()
    # nn.nn1()
    # # nn.nn1_single_edge_length(epochs=3000)
    # print(f"Time elapsed: {round(time.time() - start_time, 5)} seconds")

    # New procrustes test
    # import math

    # edge_count = 8
    # test = pp.create_random_ngon(edge_count)
    # pro = pp.procrustes(test)['transformed_contour']
    # pro_simp = pp.simple_procrustes(test)

    # pp.plot_polygon(pro, '--')
    # pp.plot_polygon(pro_simp['contour'], 'g')
    # pp.plot_polygon(pro_simp["contour"] @ pro_simp["rotation"]
    #                 * pro_simp["scale"] + pro_simp["translation"], 'r:')
    # pp.plot_polygon(pp.create_regular_ngon(edge_count), 'k')

    # plt.gca().set_aspect("equal", adjustable="box")
    # plt.show()
    pass
