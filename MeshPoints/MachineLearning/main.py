import contour_transform as ct
import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    test_contour = np.array([[2.02, 0, -2.14, -2.34, -0.24, 1.78], [1.52, 3, 1.64, -0.78, -2.18, -0.5]])

    transformed_contour = ct.procrustes_superimposition(test_contour)
    # ct.plot_ngon(test_contour)
    ct.plot_ngon(transformed_contour['transformed_contour'])
    plt.show()

    print(transformed_contour)