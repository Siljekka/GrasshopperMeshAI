import numpy as np
import csv
import pre_processing as pp
import gmsh


def generate_dataset(dataset_size: int, num_sides: int, target_edge_length: float) -> list:

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


def generate_dataset_all_lc(dataset_size: int, num_sides: int) -> list:

    # Start a gmsh API session
    gmsh.initialize()

    dataset = []
    for _ in range(dataset_size):
        test_polygon = pp.create_random_ngon(num_sides)
        procrustes = pp.procrustes(test_polygon)

        dataset_nested = pp.mesh_contour_all_lc(procrustes["transformed_contour"])
        for l in dataset_nested:
            dataset.append(l)

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

def single_edge_length_mesh_to_csv(features, dataset_size: int, number_of_sides: int) -> None:
    # Columns needed:
    #  contour nodes ( xi | yi) | target_edge_length | num inner nodes
    with open(f"data/{number_of_sides}-gon-lc-04-mesh-dataset.csv", "w", newline="") as file:
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

def mesh_to_csv(features, dataset_size: int, number_of_sides: int) -> None:
    # Columns needed:
    #  contour nodes ( xi | yi) | target_edge_length | num inner nodes
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
