import contour_transform as ct
import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    test_contour = np.array(
        [[7.06, 9.16, 9.62, 14.2, 14.36, 10.46], [0.55, -0.99, -3.25, -0.95, 3.31, 3.53]])


    transformed_contour = ct.procrustes_superimposition_scipy(test_contour)
    # ct.plot_ngon(test_contour)
    ct.plot_ngon(np.array([[x[0] for x in transformed_contour[1]], [y[1] for y in transformed_contour[1]]]))
    ct.plot_ngon(np.array([[x[0] for x in transformed_contour[0]], [y[1] for y in transformed_contour[0]]]))
    # ct.plot_ngon(transformed_contour['transformed_contour'])
    plt.show()

    print(transformed_contour)
