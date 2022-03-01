import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage as ndimage
from PIL import Image, ImageOps


class LatGeo:
    """
    LAT Geometry.
    """

    F_2 = 7000
    n_si = 3.416
    n_vac = 1.0

    th_fwhp = 35 * np.pi / 180
    n_scan = 100

    de_ang = 0.5 / 60 * np.pi / 180  # arcsec = 1/60 degree
    lambda_ = (30.0 / 150.0) * 0.01  # [m]
    k = 2 * np.pi / lambda_

    lyot_y = 453.62

    x_ap = 0
    y_ap = 0
    z_ap = 0

    y_source = 0


th1_l1 = np.pi / 2
th2_l1 = np.pi / 2
th1_l2 = np.pi / 2
th2_l2 = np.pi / 2
th1_l3 = np.pi / 2
th2_l3 = np.pi / 2

lens_t1 = 4.5633e1
lens_t2 = 4.9101e1
lens_t3 = 3.1790e1

lens3_y = (1.5 + 0.45 + 2.1553254002897502) * 1e1
lens2_y = lens3_y + 10.676699903527055e1 + lens_t3
lens1_y = lens2_y + 3.8933e1 + +0.45e1 + 0.5e1 + \
    0.45e1 + 47.776194563951336e1 + lens_t2
y_lyot = lens1_y + 1.0553506281444689e1 + lens_t1


def z1b(x, y):
    """
    Surface shape of lens 1 side B.
    """

    x_temp = x / 1e1
    y_temp = y / 1e1

    c = -5.97879325e-3
    k = 29.99835006
    a_1 = 6.86096781e-3
    a_2 = 4.44719534e-7
    a_3 = -2.62178001e-9

    r = np.sqrt(x_temp ** 2 + y_temp ** 2)
    if 1 == 0:
        amp = 0
    else:
        amp = (c * r ** 2) / (1 + np.sqrt(1 - ((1 + k) * c ** 2 * r ** 2)))
        amp += a_1 * r ** 2 + a_2 * r ** 4 + a_3 * r ** 6
    return amp * 10


def z1a(x, y):
    """
    Surface shape of lens 1 side A.
    """
    x_temp = x / 1e1
    y_temp = y / 1e1

    c = 9.65674282e-3
    k = -30.00122787
    a_1 = 1.45021174e-3
    a_2 = 2.84858123e-6
    a_3 = -5.02176737e-9

    r = np.sqrt(x_temp ** 2 + y_temp ** 2)
    if 1 == 0:
        amp = 0
    else:
        amp = (c * r ** 2) / (1 + np.sqrt(1 - ((1 + k) * c ** 2 * r ** 2)))
        amp += a_1 * r ** 2 + a_2 * r ** 4 + a_3 * r ** 6
    return amp * 10


# fix units in all the d_ functions!!!!!


def d_z1b(x, y):
    """
    Normal vector of surface of lens 1 side A.
    """
    c = -5.97879325e-3
    k = 29.99835006
    a_1 = 6.86096781e-3
    a_2 = 4.44719534e-7
    a_3 = -2.62178001e-9

    x_temp = x / 1e1
    y_temp = y / 1e1
    r = np.sqrt(x_temp ** 2 + y_temp ** 2)

    coeff_1 = (a_1 * 2) + (a_2 * 4 * r ** 2) + (a_3 * 6 * r ** 4)
    coeff_2 = (c * 2) / (1 + np.sqrt(1 - ((1 + k) * c ** 2 * r ** 2)))
    coeff_3 = (c ** 3 * (k + 1) * r ** 2) / (
        np.sqrt(1 - ((1 + k) * c ** 2 * r ** 2))
        * (1 + np.sqrt(1 - ((1 + k) * c ** 2 * r ** 2))) ** 2
    )
    amp_x = x_temp * (coeff_1 + coeff_2 + coeff_3)
    amp_y = y_temp * (coeff_1 + coeff_2 + coeff_3)

    return amp_x, amp_y


def d_z1a(x, y):
    """
    Normal vector of surface of lens 1 side A.
    """
    c = 9.65674282e-3
    k = -30.00122787
    a_1 = 1.45021174e-3
    a_2 = 2.84858123e-6
    a_3 = -5.02176737e-9

    x_temp = x / 1e1
    y_temp = y / 1e1
    r = np.sqrt(x_temp ** 2 + y_temp ** 2)

    coeff_1 = (a_1 * 2) + (a_2 * 4 * r ** 2) + (a_3 * 6 * r ** 4)
    coeff_2 = (c * 2) / (1 + np.sqrt(1 - ((1 + k) * c ** 2 * r ** 2)))
    coeff_3 = (c ** 3 * (k + 1) * r ** 2) / (
        np.sqrt(1 - ((1 + k) * c ** 2 * r ** 2))
        * (1 + np.sqrt(1 - ((1 + k) * c ** 2 * r ** 2))) ** 2
    )
    amp_x = x_temp * (coeff_1 + coeff_2 + coeff_3)
    amp_y = y_temp * (coeff_1 + coeff_2 + coeff_3)

    return amp_x, amp_y


def m1a_into_tele(x, y, z):

    x_rot2 = x
    y_rot2 = (y * np.cos(th1_l1) - z * np.sin(th1_l1)) + (lens1_y + lens_t1)
    z_rot2 = y * np.sin(th1_l1) + z * np.cos(th1_l1)
    return x_rot2, y_rot2, z_rot2


def m1b_into_tele(x, y, z):

    x_rot1 = x
    y_rot1 = y * np.cos(th1_l1) - z * np.sin(th1_l1) + lens1_y
    z_rot1 = y * np.sin(th1_l1) + z * np.cos(th1_l1)
    return x_rot1, y_rot1, z_rot1


def tele_into_m1b(x, y, z):

    y -= lens1_y
    x_temp = x
    y_temp = y * np.cos(-th1_l1) - z * np.sin(-th1_l1)
    z_temp = y * np.sin(-th1_l1) + z * np.cos(-th1_l1)

    return x_temp, y_temp, z_temp


def tele_into_m1a(x, y, z):

    y -= lens1_y + lens_t1
    x_temp = x
    y_temp = y * np.cos(-th1_l1) - z * np.sin(-th1_l1)
    z_temp = y * np.sin(-th1_l1) + z * np.cos(-th1_l1)

    return x_temp, y_temp, z_temp


def z2a(x, y):
    """
    Surface shape of lens 2 side A.
    """
    x_temp = x / 1e1
    y_temp = y / 1e1

    r = np.sqrt(x_temp ** 2 + y_temp ** 2)

    c = 1.22111813e-02
    k = -30.00046667
    a_1 = 1.56460553e-03
    a_2 = 3.06349779e-06
    a_3 = -5.00246955e-09

    if 1 == 0:
        amp = 0
    else:
        amp = (c * r ** 2) / (1 + np.sqrt(1 - (1 + k) * c ** 2 * r ** 2))
        amp += a_1 * r ** 2 + a_2 * r ** 4 + a_3 * r ** 6

    return amp * 10


def z2b(x, y):
    """
    Surface shape of lens 2 side B.
    """
    x_temp = x / 1e1
    y_temp = y / 1e1

    r = np.sqrt(x_temp ** 2 + y_temp ** 2)

    c = 2.39616010e-02
    k = -1.04603079
    a_1 = -8.38235791e-03
    a_2 = -1.39048036e-06
    a_3 = 1.84669164e-10

    if 1 == 0:
        amp = 0
    else:
        amp = (c * r ** 2) / (1 + np.sqrt(1 - (1 + k) * c ** 2 * r ** 2))
        amp += a_1 * r ** 2 + a_2 * r ** 4 + a_3 * r ** 6

    return amp * 10


def d_z2a(x, y):

    c = 1.22111813e-02
    k = -30.00046667
    a_1 = 1.56460553e-03
    a_2 = 3.06349779e-06
    a_3 = -5.00246955e-09

    x_temp = x / 1e1
    y_temp = y / 1e1

    r = np.sqrt(x_temp ** 2 + y_temp ** 2)

    coeff_1 = (a_1 * 2) + (a_2 * 4 * r ** 2) + (a_3 * 6 * r ** 4)
    coeff_2 = (c * 2) / (1 + np.sqrt(1 - ((1 + k) * c ** 2 * r ** 2)))
    coeff_3 = (c ** 3 * (k + 1) * r ** 2) / (
        np.sqrt(1 - ((1 + k) * c ** 2 * r ** 2))
        * (1 + np.sqrt(1 - ((1 + k) * c ** 2 * r ** 2))) ** 2
    )
    amp_x = x_temp * (coeff_1 + coeff_2 + coeff_3)
    amp_y = y_temp * (coeff_1 + coeff_2 + coeff_3)
    return amp_x, amp_y


def d_z2b(x, y):

    c = 2.39616010e-02
    k = -1.04603079
    a_1 = -8.38235791e-03
    a_2 = -1.39048036e-06
    a_3 = 1.84669164e-10

    x_temp = x / 1e1
    y_temp = y / 1e1

    r = np.sqrt(x_temp ** 2 + y_temp ** 2)

    coeff_1 = (a_1 * 2) + (a_2 * 4 * r ** 2) + (a_3 * 6 * r ** 4)
    coeff_2 = (c * 2) / (1 + np.sqrt(1 - ((1 + k) * c ** 2 * r ** 2)))
    coeff_3 = (c ** 3 * (k + 1) * r ** 2) / (
        np.sqrt(1 - ((1 + k) * c ** 2 * r ** 2))
        * (1 + np.sqrt(1 - ((1 + k) * c ** 2 * r ** 2))) ** 2
    )
    amp_x = x_temp * (coeff_1 + coeff_2 + coeff_3)
    amp_y = y_temp * (coeff_1 + coeff_2 + coeff_3)

    return amp_x, amp_y


def m2a_into_tele(x, y, z):

    x_rot2 = x
    y_rot2 = (y * np.cos(th2_l2) - z * np.sin(th2_l2)) + (lens2_y + 4.9101e1)
    z_rot2 = y * np.sin(th2_l2) + z * np.cos(th2_l2)
    return x_rot2, y_rot2, z_rot2


def m2b_into_tele(x, y, z):

    x_rot1 = x
    y_rot1 = y * np.cos(th1_l2) - z * np.sin(th1_l2) + lens2_y
    z_rot1 = y * np.sin(th1_l2) + z * np.cos(th1_l2)
    return x_rot1, y_rot1, z_rot1


def tele_into_m2b(x, y, z):
    y -= lens2_y
    x_temp = x
    y_temp = y * np.cos(-th1_l2) - z * np.sin(-th1_l2)
    z_temp = y * np.sin(-th1_l2) + z * np.cos(-th1_l2)

    return x_temp, y_temp, z_temp


def tele_into_m2a(x, y, z):
    y -= lens2_y + 4.9101e1
    x_temp = x
    y_temp = y * np.cos(-th2_l2) - z * np.sin(-th2_l2)
    z_temp = y * np.sin(-th2_l2) + z * np.cos(-th2_l2)

    return x_temp, y_temp, z_temp


def z3a(x, y):
    """
    Surface shape of lens 3 side A.
    """
    x_temp = x / 1e1
    y_temp = y / 1e1

    c = -2.48540190e-02
    k = -30.00640716
    a_1 = 4.05886709e-04
    a_2 = 5.00081009e-06
    a_3 = 5.00669812e-09

    r = np.sqrt(x_temp ** 2 + y_temp ** 2)
    if 1 == 0:
        amp = 0
    else:
        amp = (c * r ** 2) / (1 + np.sqrt(1 - ((1 + k) * c ** 2 * r ** 2)))
        amp += a_1 * r ** 2 + a_2 * r ** 4 + a_3 * r ** 6
    return amp * 10


def z3b(x, y):
    """
    Surface shape of lens 3 side B.
    """
    x_temp = x / 1e1
    y_temp = y / 1e1

    c = -2.50138371e-02
    k = -18.18966463
    a_1 = 2.29039839e-03
    a_2 = 5.00824363e-06
    a_3 = 1.94511404e-09

    r = np.sqrt(x_temp ** 2 + y_temp ** 2)
    if 1 == 0:
        amp = 0
    else:
        amp = (c * r ** 2) / (1 + np.sqrt(1 - ((1 + k) * c ** 2 * r ** 2)))
        amp += a_1 * r ** 2 + a_2 * r ** 4 + a_3 * r ** 6
    return amp * 10


def filt(x, y):
    """
    Surface shape of filter.
    """
    return 0


def d_z3a(x, y):
    """
    Normal vector on surface of lens 3 side A.
    """
    c = -2.48540190e-02
    k = -30.00640716
    a_1 = 4.05886709e-04
    a_2 = 5.00081009e-06
    a_3 = 5.00669812e-09

    x_temp = x / 1e1
    y_temp = y / 1e1
    r = np.sqrt(x_temp ** 2 + y_temp ** 2)

    coeff_1 = (a_1 * 2) + (a_2 * 4 * r ** 2) + (a_3 * 6 * r ** 4)
    coeff_2 = (c * 2) / (1 + np.sqrt(1 - ((1 + k) * c ** 2 * r ** 2)))
    coeff_3 = (c ** 3 * (k + 1) * r ** 2) / (
        np.sqrt(1 - ((1 + k) * c ** 2 * r ** 2))
        * (1 + np.sqrt(1 - ((1 + k) * c ** 2 * r ** 2))) ** 2
    )
    amp_x = x_temp * (coeff_1 + coeff_2 + coeff_3)
    amp_y = y_temp * (coeff_1 + coeff_2 + coeff_3)

    return amp_x, amp_y


def d_z3b(x, y):
    """
    Normal vector on surface of lens 3 side B.
    """
    x_temp = x / 1e1
    y_temp = y / 1e1

    c = -2.50138371e-02
    k = -18.18966463
    a_1 = 2.29039839e-03
    a_2 = 5.00824363e-06
    a_3 = 1.94511404e-09

    r = np.sqrt(x_temp ** 2 + y_temp ** 2)

    coeff_1 = (a_1 * 2) + (a_2 * 4 * r ** 2) + (a_3 * 6 * r ** 4)
    coeff_2 = (c * 2) / (1 + np.sqrt(1 - ((1 + k) * c ** 2 * r ** 2)))
    coeff_3 = (c ** 3 * (k + 1) * r ** 2) / (
        np.sqrt(1 - ((1 + k) * c ** 2 * r ** 2))
        * (1 + np.sqrt(1 - ((1 + k) * c ** 2 * r ** 2))) ** 2
    )
    amp_x = x_temp * (coeff_1 + coeff_2 + coeff_3)
    amp_y = y_temp * (coeff_1 + coeff_2 + coeff_3)

    return amp_x, amp_y


def m3a_into_tele(x, y, z):
    """
    Coordinate transformation from lens 3 side A to telescope.
    """
    x_rot1 = x
    y_rot1 = y * np.cos(th1_l3) - z * np.sin(th1_l3) + (lens3_y + 3.1790e1)
    z_rot1 = y * np.sin(th1_l3) + z * np.cos(th1_l3)
    return x_rot1, y_rot1, z_rot1


def m3b_into_tele(x, y, z):
    """
    Coordinate transformation from lens 3 side B to telescope.
    """
    x_rot2 = x
    y_rot2 = (y * np.cos(th2_l3) - z * np.sin(th2_l3)) + (lens3_y)
    z_rot2 = y * np.sin(th2_l3) + z * np.cos(th2_l3)
    return x_rot2, y_rot2, z_rot2


def tele_into_m3a(x, y, z):
    """
    Coordinate transformation from telescope to lens 3, side A.
    """
    y -= lens3_y + 3.1790e1
    x_temp = x
    y_temp = y * np.cos(-th1_l3) - z * np.sin(-th1_l3)
    z_temp = y * np.sin(-th1_l3) + z * np.cos(-th1_l3)

    return x_temp, y_temp, z_temp


def tele_into_m3b(x, y, z):
    """
    Coordinate transformation from telescope to lens 3, side B.
    """
    y -= lens3_y
    x_temp = x
    y_temp = y * np.cos(-th2_l3) - z * np.sin(-th2_l3)
    z_temp = y * np.sin(-th2_l3) + z * np.cos(-th2_l3)

    return x_temp, y_temp, z_temp


def plot_lenses():
    """
    Plots geometry of LATR OT lenses.
    """
    x = np.linspace(-(448 / 2), (448 / 2), 100)  # [mm]
    y = np.linspace(-(448 / 2), (448 / 2), 100)  # [mm]
    X, Y = np.meshgrid(x, y)

    r = np.sqrt(X ** 2 + Y ** 2)

    Z1a = z1a(X, Y)
    Z1b = z1b(X, Y)
    Z1a = np.where(r < 300, Z1a, np.nan)
    Z1b = np.where(r < 300, Z1b, np.nan)

    Z2a = z2a(X, Y)
    Z2b = z2b(X, Y)
    Z2a = np.where(r < 300, Z2a, np.nan)
    Z2b = np.where(r < 300, Z2b, np.nan)

    Z3a = z3a(X, Y)
    Z3b = z3b(X, Y)
    Z3a = np.where(r < 300, Z3a, np.nan)
    Z3b = np.where(r < 300, Z3b, np.nan)

    Xt1b, Yt1b, Zt1b = m1b_into_tele(X, Y, Z1b)
    Xt1a, Yt1a, Zt1a = m1a_into_tele(X, Y, Z1a)

    Xt2b, Yt2b, Zt2b = m2b_into_tele(X, Y, Z2b)
    Xt2a, Yt2a, Zt2a = m2a_into_tele(X, Y, Z2a)

    Xt3b, Yt3b, Zt3b = m3b_into_tele(X, Y, Z3b)
    Xt3a, Yt3a, Zt3a = m3a_into_tele(X, Y, Z3a)

    plt.plot(
        Yt1b[:, int(len(Yt1a) / 2)],
        Zt1b[:, int(len(Yt1a) / 2)],
        "-",
        color="k",
        label="1b",
    )
    plt.plot(
        Yt1a[:, int(len(Yt1a) / 2)],
        Zt1a[:, int(len(Yt1a) / 2)],
        "-",
        color="k",
        label="1a",
    )
    plt.plot(
        Yt2b[:, int(len(Yt2a) / 2)],
        Zt2b[:, int(len(Yt2a) / 2)],
        "-",
        color="k",
        label="2b",
    )
    plt.plot(
        Yt2a[:, int(len(Yt2a) / 2)],
        Zt2a[:, int(len(Yt2a) / 2)],
        "-",
        color="k",
        label="2a",
    )
    plt.plot(
        Yt3b[:, int(len(Yt3a) / 2)],
        Zt3b[:, int(len(Yt3a) / 2)],
        "-",
        color="k",
        label="3b",
    )

    plt.plot(
        Yt3a[:, int(len(Yt3a) / 2)],
        Zt3a[:, int(len(Yt3a) / 2)],
        "-",
        color="k",
        label="3a",
    )
