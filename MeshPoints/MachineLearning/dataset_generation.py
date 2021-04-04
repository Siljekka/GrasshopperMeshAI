import numpy as np
import csv


def contour_to_csv(triangles: np.array, number_of_sides: int) -> None:
    with open(f"data/{number_of_sides}-gon-dataset", "w", newline="") as file:
        writer = csv.writer(file)
        header = ["id"]
        for i in range(1, number_of_sides + 1):
            header.append(f"x{i}")
            header.append(f"y{i}")
        header.append("scaling_ratio")

        writer.writerow(header)


def mesh_to_csv() -> None:
    # Columns needed:
    #  edge_nodes ( xi | yi) | inner_nodes (xi | yi) | num_inner_nodes | scaling factor
    pass
