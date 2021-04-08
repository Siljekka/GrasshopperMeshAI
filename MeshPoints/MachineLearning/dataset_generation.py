import numpy as np
import csv
import pre_processing as pp
import gmsh


def generate_dataset(dataset_size: int, num_sides: int, target_edge_length: float) -> list:
    # dataset_size = {
    #     2: 3_000,
    #     3: 4_500,
    #     4: 6_000,
    #     5: 9_000,
    #     6: 12_000,
    #     7: 17_000,
    #     8: 24_000,
    #     9: 35_000,
    #     10: 48_000,
    #     12: 95_000,
    #     14: 190_000,
    #     16: 380_000
    # }

    # Start a gmsh API session
    gmsh.initialize()
    dataset = []
    for _ in range(dataset_size):
        test_polygon = pp.create_random_ngon(num_sides)
        procrustes = pp.procrustes(test_polygon)

        dataset.append(pp.mesh_contour(
            procrustes["transformed_contour"], target_edge_length))

    # Close API call
    gmsh.finalize()

    return dataset


def contour_to_csv(triangles: np.array, number_of_sides: int) -> None:
    with open(f"data/{number_of_sides}-gon-dataset", "w", newline="") as file:
        writer = csv.writer(file)
        header = ["id"]
        for i in range(1, number_of_sides + 1):
            header.append(f"x{i}")
            header.append(f"y{i}")
        header.append("scaling_ratio")

        writer.writerow(header)


def mesh_to_csv(features, dataset_size: int, number_of_sides: int) -> None:
    # Columns needed:
    #  contour nodes ( xi | yi) | target_edge_length |num inner nodes
    with open(f"data/{number_of_sides}-gon-mesh-dataset.csv", "w", newline="") as file:
        writer = csv.writer(file)
        header = []
        for i in range(1, number_of_sides + 1):
            header.append(f"x{i}")
            header.append(f"y{i}")
        header.append("target_edge_length")
        header.append("internal_nodes")

        writer.writerow(header)

        for i, contour in enumerate(features):
            print("writing", i, "of", dataset_size, end="\r")
            writer.writerow(contour)
