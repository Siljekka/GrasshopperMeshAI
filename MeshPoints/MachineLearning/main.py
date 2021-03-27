import contour_transform as ct
import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    # ngon = ct.create_regular_ngon(10)
    # ct.plot_ngon(ngon)
    test_contour = np.array([[0.82, 4.56, 3.62, -2.58, -4.28, -0.72], [-3.69, -0.45, 3.59, 2.33, -0.91, -1.25]])

    transformed_contour = ct.procrustes_superimposition(test_contour)

    ct.plot_ngon(transformed_contour)
    plt.show()