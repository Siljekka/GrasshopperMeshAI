import contour_transform as ct


if __name__ == '__main__':
    ngon = ct.create_regular_ngon(10)
    ct.plot_2d_points(ngon)
