from math import cos, sin, pi
import matplotlib.pyplot as plt


def create_regular_ngon(number_of_sides: int):
    """
    Creates a regular n-gon with circumradius = 1.
    This is equivalent to creating a regular ngon inscribed in a unit circle.
    """
    x = []
    y = []
    for n in range(number_of_sides):
        x.append(cos(2*pi*n/number_of_sides))
        y.append(sin(2*pi*n/number_of_sides))

    return [x, y]


def plot_2d_points(coords):
    plt.plot(coords[0], coords[1], 'o')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()
