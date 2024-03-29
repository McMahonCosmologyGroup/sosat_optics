import matplotlib
import matplotlib.pyplot as plt
import numpy as np
# from ot_geo import *
from scipy import interpolate, optimize
from tqdm import tqdm

from sosat_optics import ot_geo
from sosat_optics.ot_geo import *

matplotlib.rcParams["font.size"] = 16
matplotlib.rcParams.update(
    {
        "axes.grid": False,
        "grid.color": "grey",
        "grid.alpha": 0.2,
        "xtick.direction": "in",
        "ytick.direction": "in",
    }
)

red = "#e42536"
blue = "#570efc"
orange = "#f89c20"


def snell_vec(n1, n2, N_surf, s1):
    # s1 is the incoming vector, pointing from the light source to the surface
    # N_surf is the normal of the surface

    s2 = (n1 / n2) * np.cross(N_surf, (np.cross(-N_surf, s1))) - N_surf * np.sqrt(
        1 - (n1 / n2) ** 2 * np.dot((np.cross(N_surf, s1)), (np.cross(N_surf, s1)))
    )

    return s2


def aperature_fields(P_rx, tele_geo, plot, col):
    alph = 0.05  # transparency of plotted lines

    y_ap = tele_geo.y_ap
    horn_fwhp = tele_geo.th_fwhp
    n_vac = tele_geo.n_vac
    n_si = tele_geo.n_si

    N_linear = tele_geo.n_scan

    # Step 1:  grid the plane of rays shooting out of receiver feed
    theta = np.linspace((np.pi / 2) - 0.35, (np.pi / 2) + 0.35, N_linear)
    phi = np.linspace((np.pi / 2) - 0.35, (np.pi / 2) + 0.35, N_linear)

    theta, phi = np.meshgrid(theta, phi)
    theta = np.ravel(theta)
    phi = np.ravel(phi)

    # Step 2: calculate the position + local surface normal for the dish
    n_pts = len(theta)
    out = np.zeros((17, n_pts))

    for ii in range(n_pts):

        th = theta[ii]
        ph = phi[ii]

        r_hat = [np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)]

        alpha = r_hat[0]
        beta = r_hat[1]
        gamma = r_hat[2]

        # Receiver feed position [mm] (in telescope reference frame):
        x_0 = P_rx[0]
        y_0 = P_rx[1]
        z_0 = P_rx[2]

        def root_z3b(t):

            x = x_0 + alpha * t
            y = y_0 + beta * t
            z = z_0 + gamma * t

            xm3b, ym3b, zm3b = tele_into_m3b(
                x, y, z
            )  # Convert ray's endpoint into M2 coordinates

            z_m3b = z3b(xm3b, ym3b)  # Z of mirror in M2 coordinates
            if np.isnan(z_m3b):
                z_m3b = 0
            root = zm3b - z_m3b
            return root

        t_m3b = optimize.brentq(root_z3b, 2, 1600)

        # Location of where ray hits M2
        x_m3b = x_0 + alpha * t_m3b
        y_m3b = y_0 + beta * t_m3b
        z_m3b = z_0 + gamma * t_m3b
        P_m3b = np.array([x_m3b, y_m3b, z_m3b])

        ###### in M2 coordinates ##########################
        x_m3b_temp, y_m3b_temp, z_m3b_temp = tele_into_m3b(
            x_m3b, y_m3b, z_m3b
        )  # P_m2 temp
        x_rx_temp, y_rx_temp, z_rx_temp = tele_into_m3b(x_0, y_0, z_0)  # P_rx temp
        norm = d_z3b(x_m3b_temp, y_m3b_temp)
        norm_temp = np.array([-norm[0], -norm[1], 1])

        N_hat = norm_temp / np.sqrt(sum(norm_temp ** 2))

        vec_rx_m3b = np.array([x_m3b_temp, y_m3b_temp, z_m3b_temp]) - np.array(
            [x_rx_temp, y_rx_temp, z_rx_temp]
        )
        dist_rx_m3b = np.sqrt(np.sum(vec_rx_m3b ** 2))
        tan_rx_m3b = vec_rx_m3b / dist_rx_m3b

        # Use Snell's Law to find angle of outgoing ray:

        tan_og_si = snell_vec(n_vac, n_si, N_hat, tan_rx_m3b)

        # Transform back to telescope cordinates ############

        N_hat_t = np.zeros(3)
        N_hat_t[0] = N_hat[0]
        N_hat_t[1] = N_hat[1] * np.cos(th2_l3) - N_hat[2] * np.sin(th2_l3)
        N_hat_t[2] = N_hat[1] * np.sin(th2_l3) + N_hat[2] * np.cos(th2_l3)

        tan_rx_m3b_t = np.zeros(3)
        tan_rx_m3b_t[0] = tan_rx_m3b[0]
        tan_rx_m3b_t[1] = tan_rx_m3b[1] * np.cos(th2_l3) - tan_rx_m3b[2] * np.sin(
            th2_l3
        )
        tan_rx_m3b_t[2] = tan_rx_m3b[1] * np.sin(th2_l3) + tan_rx_m3b[2] * np.cos(
            th2_l3
        )

        tan_og_t = np.zeros(3)
        tan_og_t[0] = tan_og_si[0]
        tan_og_t[1] = tan_og_si[1] * np.cos(th2_l3) - tan_og_si[2] * np.sin(th2_l3)
        tan_og_t[2] = tan_og_si[1] * np.sin(th2_l3) + tan_og_si[2] * np.cos(th2_l3)

        alpha = tan_og_t[0]
        beta = tan_og_t[1]
        gamma = tan_og_t[2]

        def root_z3a(t):
            x = P_m3b[0] + alpha * t
            y = P_m3b[1] + beta * t
            z = P_m3b[2] + gamma * t
            xm3a, ym3a, zm3a = tele_into_m3a(
                x, y, z
            )  # take ray end coordinates and convert to M1 coordinates
            z_m3a = z3a(xm3a, ym3a)  # Z of mirror 1 in M1 coordinates
            if np.isnan(z_m3a):
                z_m3a = 0
            root = zm3a - z_m3a
            return root

        t_m3a = optimize.brentq(root_z3a, 5, 200)

        # Location of where ray hits M1
        x_m3a = P_m3b[0] + alpha * t_m3a
        y_m3a = P_m3b[1] + beta * t_m3a
        z_m3a = P_m3b[2] + gamma * t_m3a
        P_m3a = np.array([x_m3a, y_m3a, z_m3a])

        ###### in M1 cordinates ##########################
        x_m3a_temp, y_m3a_temp, z_m3a_temp = tele_into_m3a(x_m3a, y_m3a, z_m3a)
        x_m3b_temp, y_m3b_temp, z_m3b_temp = tele_into_m3a(
            P_m3b[0], P_m3b[1], P_m3b[2]
        )  # P_1a temp
        norm = d_z3a(x_m3a_temp, y_m3a_temp)
        norm_temp = np.array([-norm[0], -norm[1], 1])
        N_hat = norm_temp / np.sqrt(sum(norm_temp ** 2))
        vec_m3b_m3a = np.array([x_m3a_temp, y_m3a_temp, z_m3a_temp]) - np.array(
            [x_m3b_temp, y_m3b_temp, z_m3b_temp]
        )
        dist_m3b_m3a = np.sqrt(np.sum(vec_m3b_m3a ** 2))
        tan_m3b_m3a = vec_m3b_m3a / dist_m3b_m3a

        tan_og_vac = snell_vec(n_si, n_vac, N_hat, tan_m3b_m3a)

        # Transform back to telescope cordinates
        N_hat_t = np.zeros(3)
        tan_m3b_m3a_t = np.zeros(3)
        tan_og_t = np.zeros(3)

        N_hat_t[0] = N_hat[0]
        N_hat_t[1] = N_hat[1] * np.cos(th1_l3) - N_hat[2] * np.sin(th1_l3)
        N_hat_t[2] = N_hat[1] * np.sin(th1_l3) + N_hat[2] * np.cos(th1_l3)

        tan_m3b_m3a_t[0] = tan_m3b_m3a[0]
        tan_m3b_m3a_t[1] = tan_m3b_m3a[1] * np.cos(th1_l3) - tan_m3b_m3a[2] * np.sin(
            th1_l3
        )
        tan_m3b_m3a_t[2] = tan_m3b_m3a[1] * np.sin(th1_l3) + tan_m3b_m3a[2] * np.cos(
            th1_l3
        )

        tan_og_t[0] = tan_og_vac[0]
        tan_og_t[1] = tan_og_vac[1] * np.cos(th1_l3) - tan_og_vac[2] * np.sin(th1_l3)
        tan_og_t[2] = tan_og_vac[1] * np.sin(th1_l3) + tan_og_vac[2] * np.cos(th1_l3)

        alpha = tan_og_t[0]
        beta = tan_og_t[1]
        gamma = tan_og_t[2]

        def root_z2b(t):
            x = P_m3a[0] + alpha * t
            y = P_m3a[1] + beta * t
            z = P_m3a[2] + gamma * t

            xm2b, ym2b, zm2b = tele_into_m2b(
                x, y, z
            )  # take ray end coordinates and convert to M1 coordinates
            z_m2b = z2b(xm2b, ym2b)  # Z of mirror 1 in M1 coordinates
            if np.isnan(z_m2b):
                z_m2b = 0
            root = zm2b - z_m2b
            return root

        t_m2b = optimize.brentq(root_z2b, 10, 600)

        # Location of where ray hits M1
        x_m2b = P_m3a[0] + alpha * t_m2b
        y_m2b = P_m3a[1] + beta * t_m2b
        z_m2b = P_m3a[2] + gamma * t_m2b
        P_m2b = np.array([x_m2b, y_m2b, z_m2b])

        ###### in M1 cordinates ##########################
        x_m2b_temp, y_m2b_temp, z_m2b_temp = tele_into_m2b(
            x_m2b, y_m2b, z_m2b
        )  # P_m2a temp
        x_m3a_temp, y_m3a_temp, z_m3a_temp = tele_into_m2b(
            P_m3a[0], P_m3a[1], P_m3a[2]
        )  # P_1b temp
        norm = d_z2b(x_m2b_temp, y_m2b_temp)
        norm_temp = np.array([-norm[0], -norm[1], 1])
        N_hat = norm_temp / np.sqrt(sum(norm_temp ** 2))
        vec_m3a_m2b = np.array([x_m2b_temp, y_m2b_temp, z_m2b_temp]) - np.array(
            [x_m3a_temp, y_m3a_temp, z_m3a_temp]
        )
        dist_m3a_m2b = np.sqrt(np.sum(vec_m3a_m2b ** 2))
        tan_m3a_m2b = vec_m3a_m2b / dist_m3a_m2b

        # Outgoing ray
        tan_og_si = snell_vec(n_vac, n_si, N_hat, tan_m3a_m2b)

        # Transform back to telescope cordinates
        N_hat_t = np.zeros(3)
        tan_m3a_m2b_t = np.zeros(3)
        tan_og_t = np.zeros(3)

        N_hat_t[0] = N_hat[0]
        N_hat_t[1] = N_hat[1] * np.cos(th1_l2) - N_hat[2] * np.sin(th1_l2)
        N_hat_t[2] = N_hat[1] * np.sin(th1_l2) + N_hat[2] * np.cos(th1_l2)

        tan_m3a_m2b_t[0] = tan_m3a_m2b[0]
        tan_m3a_m2b_t[1] = tan_m3a_m2b[1] * np.cos(th1_l2) - tan_m3a_m2b[2] * np.sin(
            th1_l2
        )
        tan_m3a_m2b_t[2] = tan_m3a_m2b[1] * np.sin(th1_l2) + tan_m3a_m2b[2] * np.cos(
            th1_l2
        )

        tan_og_t[0] = tan_og_si[0]
        tan_og_t[1] = tan_og_si[1] * np.cos(th1_l2) - tan_og_si[2] * np.sin(th1_l2)
        tan_og_t[2] = tan_og_si[1] * np.sin(th1_l2) + tan_og_si[2] * np.cos(th1_l2)

        alpha = tan_og_t[0]
        beta = tan_og_t[1]
        gamma = tan_og_t[2]

        def root_z2a(t):
            x = P_m2b[0] + alpha * t
            y = P_m2b[1] + beta * t
            z = P_m2b[2] + gamma * t
            xm2a, ym2a, zm2a = tele_into_m2a(
                x, y, z
            )  # take ray end coordinates and convert to M1 coordinates
            z_m2a = z2a(xm2a, ym2a)  # Z of mirror 1 in M1 coordinates
            if np.isnan(z_m2a):
                z_m2a = 0
            root = zm2a - z_m2a
            return root

        t_m2a = optimize.brentq(root_z2a, 5, 400)

        # Location of where ray hits M1
        x_m2a = P_m2b[0] + alpha * t_m2a
        y_m2a = P_m2b[1] + beta * t_m2a
        z_m2a = P_m2b[2] + gamma * t_m2a
        P_m2a = np.array([x_m2a, y_m2a, z_m2a])

        ###### in M1 cordinates ##########################
        x_m2b_temp, y_m2b_temp, z_m2b_temp = tele_into_m2a(
            x_m2b, y_m2b, z_m2b
        )  # P_m1b temp
        x_m2a_temp, y_m2a_temp, z_m2a_temp = tele_into_m2a(
            P_m2a[0], P_m2a[1], P_m2a[2]
        )  # P_1a temp
        norm = d_z2a(x_m2a_temp, y_m2a_temp)
        norm_temp = np.array([-norm[0], -norm[1], 1])
        N_hat = norm_temp / np.sqrt(sum(norm_temp ** 2))
        vec_m2b_m2a = np.array([x_m2a_temp, y_m2a_temp, z_m2a_temp]) - np.array(
            [x_m2b_temp, y_m2b_temp, z_m2b_temp]
        )
        dist_m2b_m2a = np.sqrt(np.sum(vec_m2b_m2a ** 2))
        tan_m2b_m2a = vec_m2b_m2a / dist_m2b_m2a

        # Outgoing ray

        tan_og_vac = snell_vec(n_si, n_vac, N_hat, tan_m2b_m2a)

        # Transform back to telescope cordinates
        N_hat_t = np.zeros(3)
        tan_m2b_m2a_t = np.zeros(3)
        tan_og_t = np.zeros(3)

        N_hat_t[0] = N_hat[0]
        N_hat_t[1] = N_hat[1] * np.cos(th2_l2) - N_hat[2] * np.sin(th2_l2)
        N_hat_t[2] = N_hat[1] * np.sin(th2_l2) + N_hat[2] * np.cos(th2_l2)

        tan_m2b_m2a_t[0] = tan_m2b_m2a[0]
        tan_m2b_m2a_t[1] = tan_m2b_m2a[1] * np.cos(th2_l2) - tan_m2b_m2a[2] * np.sin(
            th2_l2
        )
        tan_m2b_m2a_t[2] = tan_m2b_m2a[1] * np.sin(th2_l2) + tan_m2b_m2a[2] * np.cos(
            th2_l2
        )

        tan_og_t[0] = tan_og_vac[0]
        tan_og_t[1] = tan_og_vac[1] * np.cos(th2_l2) - tan_og_vac[2] * np.sin(th2_l2)
        tan_og_t[2] = tan_og_vac[1] * np.sin(th2_l2) + tan_og_vac[2] * np.cos(th2_l2)

        alpha = tan_og_t[0]
        beta = tan_og_t[1]
        gamma = tan_og_t[2]

        def root_z1b(t):

            x = x_m2a + alpha * t
            y = y_m2a + beta * t
            z = z_m2a + gamma * t
            xm1b, ym1b, zm1b = tele_into_m1b(
                x, y, z
            )  # Convert ray's endpoint into M2 coordinates
            z_m1b = z1b(xm1b, ym1b)  # Z of mirror in M2 coordinates
            if np.isnan(z_m1b):
                z_m1b = 0
            root = zm1b - z_m1b
            return root

        t_m1b = optimize.brentq(root_z1b, 5, 1000)

        # Location of where ray hits M2
        x_m1b = x_m2a + alpha * t_m1b
        y_m1b = y_m2a + beta * t_m1b
        z_m1b = z_m2a + gamma * t_m1b
        P_m1b = np.array([x_m1b, y_m1b, z_m1b])

        ###### in M2 coordinates ##########################
        x_m1b_temp, y_m1b_temp, z_m1b_temp = tele_into_m1b(
            x_m1b, y_m1b, z_m1b
        )  # P_m2 temp
        x_m2a_temp, y_m2a_temp, z_m2a_temp = tele_into_m1b(
            x_m2a, y_m2a, z_m2a
        )  # P_rx temp
        norm = d_z1b(x_m1b_temp, y_m1b_temp)
        norm_temp = np.array([-norm[0], -norm[1], 1])
        N_hat = norm_temp / np.sqrt(sum(norm_temp ** 2))
        vec_m2a_m1b = np.array([x_m1b_temp, y_m1b_temp, z_m1b_temp]) - np.array(
            [x_m2a_temp, y_m2a_temp, z_m2a_temp]
        )
        dist_m2a_m1b = np.sqrt(np.sum(vec_m2a_m1b ** 2))
        tan_m2a_m1b = vec_m2a_m1b / dist_m2a_m1b

        # Outgoing ray
        tan_og_si = snell_vec(n_vac, n_si, N_hat, tan_m2a_m1b)

        # Outgoing
        # Transform back to telescope cordinates
        N_hat_t = np.zeros(3)
        tan_m2a_m1b_t = np.zeros(3)
        tan_og_t = np.zeros(3)

        N_hat_t[0] = N_hat[0]
        N_hat_t[1] = N_hat[1] * np.cos(th2_l1) - N_hat[2] * np.sin(th2_l1)
        N_hat_t[2] = N_hat[1] * np.sin(th2_l1) + N_hat[2] * np.cos(th2_l1)

        tan_m2a_m1b_t[0] = tan_m2a_m1b[0]
        tan_m2a_m1b_t[1] = tan_m2a_m1b[1] * np.cos(th2_l1) - tan_m2a_m1b[2] * np.sin(
            th2_l1
        )
        tan_m2a_m1b_t[2] = tan_m2a_m1b[1] * np.sin(th2_l1) + tan_m2a_m1b[2] * np.cos(
            th2_l1
        )

        tan_og_t[0] = tan_og_si[0]
        tan_og_t[1] = tan_og_si[1] * np.cos(th2_l1) - tan_og_si[2] * np.sin(th2_l1)
        tan_og_t[2] = tan_og_si[1] * np.sin(th2_l1) + tan_og_si[2] * np.cos(th2_l1)
        ##################################################

        alpha = tan_og_t[0]
        beta = tan_og_t[1]
        gamma = tan_og_t[2]

        def root_z1a(t):
            x = P_m1b[0] + alpha * t
            y = P_m1b[1] + beta * t
            z = P_m1b[2] + gamma * t
            xm1a, ym1a, zm1a = tele_into_m1a(
                x, y, z
            )  # take ray end coordinates and convert to M1 coordinates
            z_m1a = z1a(xm1a, ym1a)  # Z of mirror 1 in M1 coordinates
            if np.isnan(z_m1a):
                z_m1a = 0
            root = zm1a - z_m1a
            return root

        t_m1a = optimize.brentq(root_z1a, 5, 200)

        # Location of where ray hits M1
        x_m1a = P_m1b[0] + alpha * t_m1a
        y_m1a = P_m1b[1] + beta * t_m1a
        z_m1a = P_m1b[2] + gamma * t_m1a
        P_m1a = np.array([x_m1a, y_m1a, z_m1a])

        ###### in M1 cordinates ##########################
        x_m1b_temp, y_m1b_temp, z_m1b_temp = tele_into_m1a(
            x_m1b, y_m1b, z_m1b
        )  # P_m1b temp
        x_m1a_temp, y_m1a_temp, z_m1a_temp = tele_into_m1a(
            P_m1a[0], P_m1a[1], P_m1a[2]
        )  # P_1a temp
        norm = d_z1a(x_m1a_temp, y_m1a_temp)
        norm_temp = np.array([-norm[0], -norm[1], 1])
        N_hat = norm_temp / np.sqrt(sum(norm_temp ** 2))
        vec_m1b_m1a = np.array([x_m1a_temp, y_m1a_temp, z_m1a_temp]) - np.array(
            [x_m1b_temp, y_m1b_temp, z_m1b_temp]
        )
        dist_m1b_m1a = np.sqrt(np.dot(vec_m1b_m1a, vec_m1b_m1a))
        tan_m1b_m1a = vec_m1b_m1a / dist_m1b_m1a

        # Outgoing ray
        tan_og_vac = snell_vec(n_si, n_vac, N_hat, tan_m1b_m1a)

        # Transform back to telescope cordinates
        N_hat_t = np.zeros(3)
        tan_m1b_m1a_t = np.zeros(3)
        tan_og_t = np.zeros(3)

        N_hat_t[0] = N_hat[0]
        N_hat_t[1] = N_hat[1] * np.cos(th1_l1) - N_hat[2] * np.sin(th1_l1)
        N_hat_t[2] = N_hat[1] * np.sin(th1_l1) + N_hat[2] * np.cos(th1_l1)

        tan_m1b_m1a_t[0] = tan_m1b_m1a[0]
        tan_m1b_m1a_t[1] = tan_m1b_m1a[1] * np.cos(th1_l1) - tan_m1b_m1a[2] * np.sin(
            th1_l1
        )
        tan_m1b_m1a_t[2] = tan_m1b_m1a[1] * np.sin(th1_l1) + tan_m1b_m1a[2] * np.cos(
            th1_l1
        )

        tan_og_t[0] = tan_og_vac[0]
        tan_og_t[1] = tan_og_vac[1] * np.cos(th1_l1) - tan_og_vac[2] * np.sin(th1_l1)
        tan_og_t[2] = tan_og_vac[1] * np.sin(th1_l1) + tan_og_vac[2] * np.cos(th1_l1)

        #         ################################################
        dist_m1a_ap = abs((y_ap - P_m1a[1]) / tan_og_t[1])
        total_path_length = (
            dist_rx_m3b
            + dist_m3b_m3a * (n_si)
            + dist_m3a_m2b
            + dist_m2b_m2a * (n_si)
            + dist_m2a_m1b
            + dist_m1b_m1a * (n_si)
            + dist_m1a_ap
        )

        pos_ap = P_m1a + dist_m1a_ap * tan_og_t

        # Estimate theta
        de_ve = np.arctan(tan_rx_m3b_t[0] / (-tan_rx_m3b_t[1]))
        de_ho = np.arctan(
            tan_rx_m3b_t[2] / np.sqrt(tan_rx_m3b_t[0] ** 2 + tan_rx_m3b_t[1] ** 2)
        )

        ################################################
        if plot == 1:
            if np.mod(ii, 51) == 0:
                if ii == 51:
                    ot_geo.plot_lenses()
                alph = 0.2
                plt.plot([y_0, y_m3b], [z_0, z_m3b], "-", color=col, alpha=alph)
                plt.plot([y_m3b, y_m3a], [z_m3b, z_m3a], "-", color=col, alpha=alph)
                #                 plt.plot([y_m3a, y_m3a + (10*tan_og_t[1])], [z_m3a, z_m3a + (10*tan_og_t[2])], "-", color='k', alpha=alph)
                plt.plot([y_m3a, y_m2b], [z_m3a, z_m2b], "-", color=col, alpha=alph)
                plt.plot([y_m2b, y_m2a], [z_m2b, z_m2a], "-", color=col, alpha=alph)
                #                 plt.plot([y_m2a,y_m2a- (50*tan_m2b_m2a_t[1])], [z_m2a,z_m2a- (50*tan_m2b_m2a_t[2])], "-", color='k')

                #                 plt.plot([y_m2a, y_m2b], [z_m2a, z_m2b], "-", color=col, alpha=alph)
                plt.plot([y_m2a, y_m1b], [z_m2a, z_m1b], "-", color=col, alpha=alph)

                plt.plot([y_m1b, y_m1a], [z_m1b, z_m1a], "-", color=col, alpha=alph)
                plt.plot(
                    [y_m1a, pos_ap[1]], [z_m1a, pos_ap[2]], "-", color=col, alpha=alph
                )

        #         Write out
        out[0, ii] = pos_ap[0]
        out[1, ii] = pos_ap[1]
        out[2, ii] = pos_ap[2]

        out[3, ii] = total_path_length
        out[4, ii] = np.exp(
            (-0.5)
            * (de_ho ** 2 + de_ve ** 2)
            / (horn_fwhp / (np.sqrt(8 * np.log(2)))) ** 2
        )

        out[5, ii] = N_hat_t[0]
        out[6, ii] = N_hat_t[1]
        out[7, ii] = N_hat_t[2]

        out[8, ii] = tan_og_t[0]
        out[9, ii] = tan_og_t[1]
        out[10, ii] = tan_og_t[2]
    return out


def rx_to_lyot_OLD(P_rx, tele_geo, plot, col):

    if plot == 1:
        ot_geo.plot_lenses()

    alph = 0.05  # transparency of plotted lines

    horn_fwhp = tele_geo.th_fwhp
    n_vac = tele_geo.n_vac
    n_si = tele_geo.n_si

    N_linear = tele_geo.n_scan

    # Step 1:  grid the plane of rays shooting out of receiver feed
    theta = np.linspace((np.pi / 2) - 0.4, (np.pi / 2) + 0.4, N_linear)
    phi = np.linspace((np.pi / 2) - 0.4, (np.pi / 2) + 0.4, N_linear)

    theta, phi = np.meshgrid(theta, phi)
    theta = np.ravel(theta)
    phi = np.ravel(phi)

    # Step 2: calculate the position + local surface normal for the dish
    n_pts = len(theta)
    out = np.zeros((17, n_pts))

    pbar = tqdm(total=n_pts)
    for ii in range(n_pts):

        th = theta[ii]
        ph = phi[ii]

        r_hat = [np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)]

        alpha = r_hat[0]
        beta = r_hat[1]
        gamma = r_hat[2]

        # Receiver feed position [mm] (in telescope reference frame):
        x_0 = P_rx[0]
        y_0 = P_rx[1]
        z_0 = P_rx[2]

        def root_z3b(t):

            x = x_0 + alpha * t
            y = y_0 + beta * t
            z = z_0 + gamma * t

            xm3b, ym3b, zm3b = tele_into_m3b(
                x, y, z
            )  # Convert ray's endpoint into M2 coordinates

            z_m3b = z3b(xm3b, ym3b)  # Z of mirror in M2 coordinates
            if np.isnan(z_m3b):
                z_m3b = 0
            root = zm3b - z_m3b
            return root

        t_m3b = optimize.brentq(root_z3b, 2, 1600)

        # Location of where ray hits M2
        x_m3b = x_0 + alpha * t_m3b
        y_m3b = y_0 + beta * t_m3b
        z_m3b = z_0 + gamma * t_m3b
        P_m3b = np.array([x_m3b, y_m3b, z_m3b])

        if x_m3b ** 2 + z_m3b ** 2 >= (392 / 2) ** 2:
            pbar.update(1)
            continue
        ###### in M2 coordinates ##########################
        x_m3b_temp, y_m3b_temp, z_m3b_temp = tele_into_m3b(
            x_m3b, y_m3b, z_m3b
        )  # P_m2 temp
        x_rx_temp, y_rx_temp, z_rx_temp = tele_into_m3b(x_0, y_0, z_0)  # P_rx temp
        norm = d_z3b(x_m3b_temp, y_m3b_temp)
        norm_temp = np.array([-norm[0], -norm[1], 1])

        N_hat = norm_temp / np.sqrt(sum(norm_temp ** 2))

        vec_rx_m3b = np.array([x_m3b_temp, y_m3b_temp, z_m3b_temp]) - np.array(
            [x_rx_temp, y_rx_temp, z_rx_temp]
        )
        dist_rx_m3b = np.sqrt(np.sum(vec_rx_m3b ** 2))
        tan_rx_m3b = vec_rx_m3b / dist_rx_m3b

        # Use Snell's Law to find angle of outgoing ray:

        tan_og_si = snell_vec(n_vac, n_si, N_hat, tan_rx_m3b)

        # Transform back to telescope cordinates ############
        N_hat_t = np.zeros(3)
        N_hat_t[0] = N_hat[0]
        N_hat_t[1] = N_hat[1] * np.cos(th2_l3) - N_hat[2] * np.sin(th2_l3)
        N_hat_t[2] = N_hat[1] * np.sin(th2_l3) + N_hat[2] * np.cos(th2_l3)

        tan_rx_m3b_t = np.zeros(3)
        tan_rx_m3b_t[0] = tan_rx_m3b[0]
        tan_rx_m3b_t[1] = tan_rx_m3b[1] * np.cos(th2_l3) - tan_rx_m3b[2] * np.sin(
            th2_l3
        )
        tan_rx_m3b_t[2] = tan_rx_m3b[1] * np.sin(th2_l3) + tan_rx_m3b[2] * np.cos(
            th2_l3
        )

        tan_og_t = np.zeros(3)
        tan_og_t[0] = tan_og_si[0]
        tan_og_t[1] = tan_og_si[1] * np.cos(th2_l3) - tan_og_si[2] * np.sin(th2_l3)
        tan_og_t[2] = tan_og_si[1] * np.sin(th2_l3) + tan_og_si[2] * np.cos(th2_l3)

        alpha = tan_og_t[0]
        beta = tan_og_t[1]
        gamma = tan_og_t[2]

        def root_z3a(t):
            x = P_m3b[0] + alpha * t
            y = P_m3b[1] + beta * t
            z = P_m3b[2] + gamma * t
            xm3a, ym3a, zm3a = tele_into_m3a(
                x, y, z
            )  # take ray end coordinates and convert to M1 coordinates
            z_m3a = z3a(xm3a, ym3a)  # Z of mirror 1 in M1 coordinates
            if np.isnan(z_m3a):
                z_m3a = 0
            root = zm3a - z_m3a
            return root

        t_m3a = optimize.brentq(root_z3a, 5, 200)

        # Location of where ray hits M1
        x_m3a = P_m3b[0] + alpha * t_m3a
        y_m3a = P_m3b[1] + beta * t_m3a
        z_m3a = P_m3b[2] + gamma * t_m3a
        P_m3a = np.array([x_m3a, y_m3a, z_m3a])

        ###### in M1 cordinates ##########################
        x_m3a_temp, y_m3a_temp, z_m3a_temp = tele_into_m3a(x_m3a, y_m3a, z_m3a)
        x_m3b_temp, y_m3b_temp, z_m3b_temp = tele_into_m3a(
            P_m3b[0], P_m3b[1], P_m3b[2]
        )  # P_1a temp
        norm = d_z3a(x_m3a_temp, y_m3a_temp)
        norm_temp = np.array([-norm[0], -norm[1], 1])
        N_hat = norm_temp / np.sqrt(sum(norm_temp ** 2))
        vec_m3b_m3a = np.array([x_m3a_temp, y_m3a_temp, z_m3a_temp]) - np.array(
            [x_m3b_temp, y_m3b_temp, z_m3b_temp]
        )
        dist_m3b_m3a = np.sqrt(np.sum(vec_m3b_m3a ** 2))
        tan_m3b_m3a = vec_m3b_m3a / dist_m3b_m3a

        tan_og_vac = snell_vec(n_si, n_vac, N_hat, tan_m3b_m3a)

        # Transform back to telescope cordinates
        N_hat_t = np.zeros(3)
        tan_m3b_m3a_t = np.zeros(3)
        tan_og_t = np.zeros(3)

        N_hat_t[0] = N_hat[0]
        N_hat_t[1] = N_hat[1] * np.cos(th1_l3) - N_hat[2] * np.sin(th1_l3)
        N_hat_t[2] = N_hat[1] * np.sin(th1_l3) + N_hat[2] * np.cos(th1_l3)

        tan_m3b_m3a_t[0] = tan_m3b_m3a[0]
        tan_m3b_m3a_t[1] = tan_m3b_m3a[1] * np.cos(th1_l3) - tan_m3b_m3a[2] * np.sin(
            th1_l3
        )
        tan_m3b_m3a_t[2] = tan_m3b_m3a[1] * np.sin(th1_l3) + tan_m3b_m3a[2] * np.cos(
            th1_l3
        )

        tan_og_t[0] = tan_og_vac[0]
        tan_og_t[1] = tan_og_vac[1] * np.cos(th1_l3) - tan_og_vac[2] * np.sin(th1_l3)
        tan_og_t[2] = tan_og_vac[1] * np.sin(th1_l3) + tan_og_vac[2] * np.cos(th1_l3)

        alpha = tan_og_t[0]
        beta = tan_og_t[1]
        gamma = tan_og_t[2]

        def root_z2b(t):
            x = P_m3a[0] + alpha * t
            y = P_m3a[1] + beta * t
            z = P_m3a[2] + gamma * t

            xm2b, ym2b, zm2b = tele_into_m2b(
                x, y, z
            )  # take ray end coordinates and convert to M1 coordinates
            z_m2b = z2b(xm2b, ym2b)  # Z of mirror 1 in M1 coordinates
            if np.isnan(z_m2b):
                z_m2b = 0
            root = zm2b - z_m2b
            return root

        t_m2b = optimize.brentq(root_z2b, 10, 600)

        # Location of where ray hits M1
        x_m2b = P_m3a[0] + alpha * t_m2b
        y_m2b = P_m3a[1] + beta * t_m2b
        z_m2b = P_m3a[2] + gamma * t_m2b
        P_m2b = np.array([x_m2b, y_m2b, z_m2b])

        ###### in M1 cordinates ##########################
        x_m2b_temp, y_m2b_temp, z_m2b_temp = tele_into_m2b(
            x_m2b, y_m2b, z_m2b
        )  # P_m2a temp
        x_m3a_temp, y_m3a_temp, z_m3a_temp = tele_into_m2b(
            P_m3a[0], P_m3a[1], P_m3a[2]
        )  # P_1b temp
        norm = d_z2b(x_m2b_temp, y_m2b_temp)
        norm_temp = np.array([-norm[0], -norm[1], 1])
        N_hat = norm_temp / np.sqrt(sum(norm_temp ** 2))
        vec_m3a_m2b = np.array([x_m2b_temp, y_m2b_temp, z_m2b_temp]) - np.array(
            [x_m3a_temp, y_m3a_temp, z_m3a_temp]
        )
        dist_m3a_m2b = np.sqrt(np.sum(vec_m3a_m2b ** 2))
        tan_m3a_m2b = vec_m3a_m2b / dist_m3a_m2b

        # Outgoing ray
        tan_og_si = snell_vec(n_vac, n_si, N_hat, tan_m3a_m2b)

        # Transform back to telescope cordinates
        N_hat_t = np.zeros(3)
        tan_m3a_m2b_t = np.zeros(3)
        tan_og_t = np.zeros(3)

        N_hat_t[0] = N_hat[0]
        N_hat_t[1] = N_hat[1] * np.cos(th1_l2) - N_hat[2] * np.sin(th1_l2)
        N_hat_t[2] = N_hat[1] * np.sin(th1_l2) + N_hat[2] * np.cos(th1_l2)

        tan_m3a_m2b_t[0] = tan_m3a_m2b[0]
        tan_m3a_m2b_t[1] = tan_m3a_m2b[1] * np.cos(th1_l2) - tan_m3a_m2b[2] * np.sin(
            th1_l2
        )
        tan_m3a_m2b_t[2] = tan_m3a_m2b[1] * np.sin(th1_l2) + tan_m3a_m2b[2] * np.cos(
            th1_l2
        )

        tan_og_t[0] = tan_og_si[0]
        tan_og_t[1] = tan_og_si[1] * np.cos(th1_l2) - tan_og_si[2] * np.sin(th1_l2)
        tan_og_t[2] = tan_og_si[1] * np.sin(th1_l2) + tan_og_si[2] * np.cos(th1_l2)

        alpha = tan_og_t[0]
        beta = tan_og_t[1]
        gamma = tan_og_t[2]

        def root_z2a(t):
            x = P_m2b[0] + alpha * t
            y = P_m2b[1] + beta * t
            z = P_m2b[2] + gamma * t
            xm2a, ym2a, zm2a = tele_into_m2a(
                x, y, z
            )  # take ray end coordinates and convert to M1 coordinates
            z_m2a = z2a(xm2a, ym2a)  # Z of mirror 1 in M1 coordinates
            if np.isnan(z_m2a):
                z_m2a = 0
            root = zm2a - z_m2a
            return root

        t_m2a = optimize.brentq(root_z2a, 5, 400)

        # Location of where ray hits M1
        x_m2a = P_m2b[0] + alpha * t_m2a
        y_m2a = P_m2b[1] + beta * t_m2a
        z_m2a = P_m2b[2] + gamma * t_m2a
        P_m2a = np.array([x_m2a, y_m2a, z_m2a])

        ###### in M1 cordinates ##########################
        x_m2b_temp, y_m2b_temp, z_m2b_temp = tele_into_m2a(
            x_m2b, y_m2b, z_m2b
        )  # P_m1b temp
        x_m2a_temp, y_m2a_temp, z_m2a_temp = tele_into_m2a(
            P_m2a[0], P_m2a[1], P_m2a[2]
        )  # P_1a temp
        norm = d_z2a(x_m2a_temp, y_m2a_temp)
        norm_temp = np.array([-norm[0], -norm[1], 1])
        N_hat = norm_temp / np.sqrt(sum(norm_temp ** 2))
        vec_m2b_m2a = np.array([x_m2a_temp, y_m2a_temp, z_m2a_temp]) - np.array(
            [x_m2b_temp, y_m2b_temp, z_m2b_temp]
        )
        dist_m2b_m2a = np.sqrt(np.sum(vec_m2b_m2a ** 2))
        tan_m2b_m2a = vec_m2b_m2a / dist_m2b_m2a

        # Outgoing ray

        tan_og_vac = snell_vec(n_si, n_vac, N_hat, tan_m2b_m2a)

        # Transform back to telescope cordinates
        N_hat_t = np.zeros(3)
        tan_m2b_m2a_t = np.zeros(3)
        tan_og_t = np.zeros(3)

        N_hat_t[0] = N_hat[0]
        N_hat_t[1] = N_hat[1] * np.cos(th2_l2) - N_hat[2] * np.sin(th2_l2)
        N_hat_t[2] = N_hat[1] * np.sin(th2_l2) + N_hat[2] * np.cos(th2_l2)

        tan_m2b_m2a_t[0] = tan_m2b_m2a[0]
        tan_m2b_m2a_t[1] = tan_m2b_m2a[1] * np.cos(th2_l2) - tan_m2b_m2a[2] * np.sin(
            th2_l2
        )
        tan_m2b_m2a_t[2] = tan_m2b_m2a[1] * np.sin(th2_l2) + tan_m2b_m2a[2] * np.cos(
            th2_l2
        )

        tan_og_t[0] = tan_og_vac[0]
        tan_og_t[1] = tan_og_vac[1] * np.cos(th2_l2) - tan_og_vac[2] * np.sin(th2_l2)
        tan_og_t[2] = tan_og_vac[1] * np.sin(th2_l2) + tan_og_vac[2] * np.cos(th2_l2)

        alpha = tan_og_t[0]
        beta = tan_og_t[1]
        gamma = tan_og_t[2]

        def root_z1b(t):

            x = x_m2a + alpha * t
            y = y_m2a + beta * t
            z = z_m2a + gamma * t
            xm1b, ym1b, zm1b = tele_into_m1b(
                x, y, z
            )  # Convert ray's endpoint into M2 coordinates
            z_m1b = z1b(xm1b, ym1b)  # Z of mirror in M2 coordinates
            if np.isnan(z_m1b):
                z_m1b = 0
            root = zm1b - z_m1b
            return root

        t_m1b = optimize.brentq(root_z1b, 5, 1000)

        # Location of where ray hits M2
        x_m1b = x_m2a + alpha * t_m1b
        y_m1b = y_m2a + beta * t_m1b
        z_m1b = z_m2a + gamma * t_m1b
        P_m1b = np.array([x_m1b, y_m1b, z_m1b])

        ###### in M2 coordinates ##########################
        x_m1b_temp, y_m1b_temp, z_m1b_temp = tele_into_m1b(
            x_m1b, y_m1b, z_m1b
        )  # P_m2 temp
        x_m2a_temp, y_m2a_temp, z_m2a_temp = tele_into_m1b(
            x_m2a, y_m2a, z_m2a
        )  # P_rx temp
        norm = d_z1b(x_m1b_temp, y_m1b_temp)
        norm_temp = np.array([-norm[0], -norm[1], 1])
        N_hat = norm_temp / np.sqrt(sum(norm_temp ** 2))
        vec_m2a_m1b = np.array([x_m1b_temp, y_m1b_temp, z_m1b_temp]) - np.array(
            [x_m2a_temp, y_m2a_temp, z_m2a_temp]
        )
        dist_m2a_m1b = np.sqrt(np.sum(vec_m2a_m1b ** 2))
        tan_m2a_m1b = vec_m2a_m1b / dist_m2a_m1b

        # Outgoing ray
        tan_og_si = snell_vec(n_vac, n_si, N_hat, tan_m2a_m1b)

        # Outgoing
        # Transform back to telescope cordinates
        N_hat_t = np.zeros(3)
        tan_m2a_m1b_t = np.zeros(3)
        tan_og_t = np.zeros(3)

        N_hat_t[0] = N_hat[0]
        N_hat_t[1] = N_hat[1] * np.cos(th2_l1) - N_hat[2] * np.sin(th2_l1)
        N_hat_t[2] = N_hat[1] * np.sin(th2_l1) + N_hat[2] * np.cos(th2_l1)

        tan_m2a_m1b_t[0] = tan_m2a_m1b[0]
        tan_m2a_m1b_t[1] = tan_m2a_m1b[1] * np.cos(th2_l1) - tan_m2a_m1b[2] * np.sin(
            th2_l1
        )
        tan_m2a_m1b_t[2] = tan_m2a_m1b[1] * np.sin(th2_l1) + tan_m2a_m1b[2] * np.cos(
            th2_l1
        )

        tan_og_t[0] = tan_og_si[0]
        tan_og_t[1] = tan_og_si[1] * np.cos(th2_l1) - tan_og_si[2] * np.sin(th2_l1)
        tan_og_t[2] = tan_og_si[1] * np.sin(th2_l1) + tan_og_si[2] * np.cos(th2_l1)
        ##################################################

        alpha = tan_og_t[0]
        beta = tan_og_t[1]
        gamma = tan_og_t[2]

        def root_z1a(t):
            x = P_m1b[0] + alpha * t
            y = P_m1b[1] + beta * t
            z = P_m1b[2] + gamma * t
            xm1a, ym1a, zm1a = tele_into_m1a(
                x, y, z
            )  # take ray end coordinates and convert to M1 coordinates
            z_m1a = z1a(xm1a, ym1a)  # Z of mirror 1 in M1 coordinates
            if np.isnan(z_m1a):
                z_m1a = 0
            root = zm1a - z_m1a
            return root

        t_m1a = optimize.brentq(root_z1a, 5, 220)

        # Location of where ray hits M1
        x_m1a = P_m1b[0] + alpha * t_m1a
        y_m1a = P_m1b[1] + beta * t_m1a
        z_m1a = P_m1b[2] + gamma * t_m1a
        P_m1a = np.array([x_m1a, y_m1a, z_m1a])

        ###### in M1 cordinates ##########################
        x_m1b_temp, y_m1b_temp, z_m1b_temp = tele_into_m1a(
            x_m1b, y_m1b, z_m1b
        )  # P_m1b temp
        x_m1a_temp, y_m1a_temp, z_m1a_temp = tele_into_m1a(
            P_m1a[0], P_m1a[1], P_m1a[2]
        )  # P_1a temp
        norm = d_z1a(x_m1a_temp, y_m1a_temp)
        norm_temp = np.array([-norm[0], -norm[1], 1])
        N_hat = norm_temp / np.sqrt(sum(norm_temp ** 2))
        vec_m1b_m1a = np.array([x_m1a_temp, y_m1a_temp, z_m1a_temp]) - np.array(
            [x_m1b_temp, y_m1b_temp, z_m1b_temp]
        )
        dist_m1b_m1a = np.sqrt(np.dot(vec_m1b_m1a, vec_m1b_m1a))
        tan_m1b_m1a = vec_m1b_m1a / dist_m1b_m1a

        # Outgoing ray
        tan_og_vac = snell_vec(n_si, n_vac, N_hat, tan_m1b_m1a)

        # Transform back to telescope cordinates
        N_hat_t = np.zeros(3)
        tan_m1b_m1a_t = np.zeros(3)
        tan_og_t = np.zeros(3)

        N_hat_t[0] = N_hat[0]
        N_hat_t[1] = N_hat[1] * np.cos(th1_l1) - N_hat[2] * np.sin(th1_l1)
        N_hat_t[2] = N_hat[1] * np.sin(th1_l1) + N_hat[2] * np.cos(th1_l1)

        tan_m1b_m1a_t[0] = tan_m1b_m1a[0]
        tan_m1b_m1a_t[1] = tan_m1b_m1a[1] * np.cos(th1_l1) - tan_m1b_m1a[2] * np.sin(
            th1_l1
        )
        tan_m1b_m1a_t[2] = tan_m1b_m1a[1] * np.sin(th1_l1) + tan_m1b_m1a[2] * np.cos(
            th1_l1
        )

        tan_og_t[0] = tan_og_vac[0]
        tan_og_t[1] = tan_og_vac[1] * np.cos(th1_l1) - tan_og_vac[2] * np.sin(th1_l1)
        tan_og_t[2] = tan_og_vac[1] * np.sin(th1_l1) + tan_og_vac[2] * np.cos(th1_l1)

        #         ################################################

        dist_m1a_lyot = abs((y_lyot - P_m1a[1]) / tan_og_t[1])
        pos_lyot = P_m1a + dist_m1a_lyot * tan_og_t

        dist_lyot_ap = abs((tele_geo.y_source - pos_lyot[1]) / tan_og_t[1])

        total_path_length = (
            dist_rx_m3b
            + dist_m3b_m3a * (n_si)
            + dist_m3a_m2b
            + dist_m2b_m2a * (n_si)
            + dist_m2a_m1b
            + dist_m1b_m1a * (n_si)
            + dist_m1a_lyot
            + dist_lyot_ap
        )

        pos_ap = pos_lyot + dist_lyot_ap * tan_og_t

        # Estimate theta
        de_ve = np.arctan(tan_rx_m3b_t[0] / (-tan_rx_m3b_t[1]))
        de_ho = np.arctan(
            tan_rx_m3b_t[2] / np.sqrt(tan_rx_m3b_t[0] ** 2 + tan_rx_m3b_t[1] ** 2)
        )

        ################################################

        if ((pos_lyot[0] ** 2 + pos_lyot[2] ** 2) <= 210 ** 2) and (
            (x_m2a ** 2 + z_m2a ** 2) <= (210 ** 2)
        ):

            if plot == 1:
                if np.mod(ii, 53) == 0:
                    alph = 0.2
                    plt.plot([y_0, y_m3b], [z_0, z_m3b], "-", color=col, alpha=alph)
                    plt.plot([y_m3b, y_m3a], [z_m3b, z_m3a], "-", color=col, alpha=alph)
                    plt.plot([y_m3a, y_m2b], [z_m3a, z_m2b], "-", color=col, alpha=alph)
                    plt.plot([y_m2b, y_m2a], [z_m2b, z_m2a], "-", color=col, alpha=alph)
                    plt.plot([y_m2a, y_m1b], [z_m2a, z_m1b], "-", color=col, alpha=alph)
                    plt.plot([y_m1b, y_m1a], [z_m1b, z_m1a], "-", color=col, alpha=alph)
                    plt.plot(
                        [y_m1a, pos_lyot[1]],
                        [z_m1a, pos_lyot[2]],
                        "-",
                        color=col,
                        alpha=alph,
                    )
                    plt.plot(
                        [pos_lyot[1], pos_ap[1]],
                        [pos_lyot[2], pos_ap[2]],
                        "-",
                        color=col,
                        alpha=alph,
                    )

        #         Write out
        out[0, ii] = pos_ap[0]
        out[1, ii] = pos_ap[1]
        out[2, ii] = pos_ap[2]

        if ((pos_lyot[0] ** 2 + pos_lyot[2] ** 2) <= (420 / 2) ** 2) and (
            (x_m2a ** 2 + z_m2a ** 2) <= (210 ** 2)
        ):
            out[3, ii] = total_path_length
            out[4, ii] = np.exp(
                (-0.5)
                * (de_ho ** 2 + de_ve ** 2)
                / (horn_fwhp / (np.sqrt(8 * np.log(2)))) ** 2
            )
        else:
            out[3, ii] = 0
            out[4, ii] = 0

        out[5, ii] = N_hat_t[0]
        out[6, ii] = N_hat_t[1]
        out[7, ii] = N_hat_t[2]

        out[8, ii] = tan_og_t[0]
        out[9, ii] = tan_og_t[1]
        out[10, ii] = tan_og_t[2]
        pbar.update(1)

    pbar.close()
    # return out

    # len_sim = int(np.sqrt(len(out[0])))
    # x_sim = np.reshape(out[0], (len_sim, len_sim))  # [mm]
    # y_sim = np.reshape(out[2], (len_sim, len_sim))  # [mm]

    # indx_x = np.where((np.isnan(x_sim) == False) & (abs(x_sim) <= 250))
    # indx_y = np.where((np.isnan(y_sim) == False) & (abs(y_sim) <= 250))

    # x_sim_new = np.linspace(np.min(x_sim[indx_x]), np.max(x_sim[indx_x]), len_sim)
    # y_sim_new = np.linspace(np.min(y_sim[indx_y]), np.max(y_sim[indx_y]), len_sim)
    # x_sim, y_sim = np.meshgrid(x_sim_new, y_sim_new)

    # a_sim = np.reshape(out[4], (len_sim, len_sim))
    # p_sim = np.reshape(out[3], (len_sim, len_sim))

    return out  # x_sim, y_sim, a_sim, p_sim


def rx_to_lyot(P_rx, tele_geo, plot, col):

    if plot:
        ot_geo.plot_lenses()

    alph = 0.05  # transparency of plotted lines

    horn_fwhp_x = tele_geo.th_fwhp_x
    horn_fwhp_y = tele_geo.th_fwhp_y

    n_vac = tele_geo.n_vac
    n_si = tele_geo.n_si

    N_linear = tele_geo.n_scan

    # Step 1:  grid the plane of rays shooting out of receiver feed
    theta = np.linspace((np.pi / 2) - 0.4, (np.pi / 2) + 0.4, N_linear)
    phi = np.linspace((np.pi / 2) - 0.4, (np.pi / 2) + 0.4, N_linear)

    theta, phi = np.meshgrid(theta, phi)
    theta = np.ravel(theta)
    phi = np.ravel(phi)

    # Step 2: calculate the position + local surface normal for the dish
    n_pts = len(theta)
    out = np.zeros((17, n_pts))

    pbar = tqdm(total=n_pts)
    for ii in range(n_pts):

        th = theta[ii]
        ph = phi[ii]

        r_hat = [np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)]

        alpha = r_hat[0]
        beta = r_hat[1]
        gamma = r_hat[2]

        # Receiver feed position [mm] (in telescope reference frame):
        x_0 = P_rx[0]
        y_0 = P_rx[1]
        z_0 = P_rx[2]

        def root_z3b(t):

            x = x_0 + alpha * t
            y = y_0 + beta * t
            z = z_0 + gamma * t

            xm3b, ym3b, zm3b = tele_into_m3b(
                x, y, z
            )  # Convert ray's endpoint into M2 coordinates

            z_m3b = z3b(xm3b, ym3b)  # Z of mirror in M2 coordinates
            if np.isnan(z_m3b):
                z_m3b = 0
            root = zm3b - z_m3b
            return root

        t_m3b = optimize.brentq(root_z3b, 2, 1600)

        # Location of where ray hits M2
        x_m3b = x_0 + alpha * t_m3b
        y_m3b = y_0 + beta * t_m3b
        z_m3b = z_0 + gamma * t_m3b
        P_m3b = np.array([x_m3b, y_m3b, z_m3b])

        if x_m3b ** 2 + z_m3b ** 2 >= (392 / 2) ** 2:
            pbar.update(1)
            continue
        ###### in M2 coordinates ##########################
        x_m3b_temp, y_m3b_temp, z_m3b_temp = tele_into_m3b(
            x_m3b, y_m3b, z_m3b
        )  # P_m2 temp
        x_rx_temp, y_rx_temp, z_rx_temp = tele_into_m3b(x_0, y_0, z_0)  # P_rx temp
        norm = d_z3b(x_m3b_temp, y_m3b_temp)
        norm_temp = np.array([-norm[0], -norm[1], 1])

        N_hat = norm_temp / np.sqrt(sum(norm_temp ** 2))

        vec_rx_m3b = np.array([x_m3b_temp, y_m3b_temp, z_m3b_temp]) - np.array(
            [x_rx_temp, y_rx_temp, z_rx_temp]
        )
        dist_rx_m3b = np.sqrt(np.sum(vec_rx_m3b ** 2))
        tan_rx_m3b = vec_rx_m3b / dist_rx_m3b

        # Use Snell's Law to find angle of outgoing ray:

        tan_og_si = snell_vec(n_vac, n_si, N_hat, tan_rx_m3b)

        # Transform back to telescope cordinates ############
        N_hat_t = np.zeros(3)
        N_hat_t[0] = N_hat[0]
        N_hat_t[1] = N_hat[1] * np.cos(th2_l3) - N_hat[2] * np.sin(th2_l3)
        N_hat_t[2] = N_hat[1] * np.sin(th2_l3) + N_hat[2] * np.cos(th2_l3)

        tan_rx_m3b_t = np.zeros(3)
        tan_rx_m3b_t[0] = tan_rx_m3b[0]
        tan_rx_m3b_t[1] = tan_rx_m3b[1] * np.cos(th2_l3) - tan_rx_m3b[2] * np.sin(
            th2_l3
        )
        tan_rx_m3b_t[2] = tan_rx_m3b[1] * np.sin(th2_l3) + tan_rx_m3b[2] * np.cos(
            th2_l3
        )

        tan_og_t = np.zeros(3)
        tan_og_t[0] = tan_og_si[0]
        tan_og_t[1] = tan_og_si[1] * np.cos(th2_l3) - tan_og_si[2] * np.sin(th2_l3)
        tan_og_t[2] = tan_og_si[1] * np.sin(th2_l3) + tan_og_si[2] * np.cos(th2_l3)

        alpha = tan_og_t[0]
        beta = tan_og_t[1]
        gamma = tan_og_t[2]

        def root_z3a(t):
            x = P_m3b[0] + alpha * t
            y = P_m3b[1] + beta * t
            z = P_m3b[2] + gamma * t
            xm3a, ym3a, zm3a = tele_into_m3a(
                x, y, z
            )  # take ray end coordinates and convert to M1 coordinates
            z_m3a = z3a(xm3a, ym3a)  # Z of mirror 1 in M1 coordinates
            if np.isnan(z_m3a):
                z_m3a = 0
            root = zm3a - z_m3a
            return root

        t_m3a = optimize.brentq(root_z3a, 5, 200)

        # Location of where ray hits M1
        x_m3a = P_m3b[0] + alpha * t_m3a
        y_m3a = P_m3b[1] + beta * t_m3a
        z_m3a = P_m3b[2] + gamma * t_m3a
        P_m3a = np.array([x_m3a, y_m3a, z_m3a])

        ###### in M1 cordinates ##########################
        x_m3a_temp, y_m3a_temp, z_m3a_temp = tele_into_m3a(x_m3a, y_m3a, z_m3a)
        x_m3b_temp, y_m3b_temp, z_m3b_temp = tele_into_m3a(
            P_m3b[0], P_m3b[1], P_m3b[2]
        )  # P_1a temp
        norm = d_z3a(x_m3a_temp, y_m3a_temp)
        norm_temp = np.array([-norm[0], -norm[1], 1])
        N_hat = norm_temp / np.sqrt(sum(norm_temp ** 2))
        vec_m3b_m3a = np.array([x_m3a_temp, y_m3a_temp, z_m3a_temp]) - np.array(
            [x_m3b_temp, y_m3b_temp, z_m3b_temp]
        )
        dist_m3b_m3a = np.sqrt(np.sum(vec_m3b_m3a ** 2))
        tan_m3b_m3a = vec_m3b_m3a / dist_m3b_m3a

        tan_og_vac = snell_vec(n_si, n_vac, N_hat, tan_m3b_m3a)

        # Transform back to telescope cordinates
        N_hat_t = np.zeros(3)
        tan_m3b_m3a_t = np.zeros(3)
        tan_og_t = np.zeros(3)

        N_hat_t[0] = N_hat[0]
        N_hat_t[1] = N_hat[1] * np.cos(th1_l3) - N_hat[2] * np.sin(th1_l3)
        N_hat_t[2] = N_hat[1] * np.sin(th1_l3) + N_hat[2] * np.cos(th1_l3)

        tan_m3b_m3a_t[0] = tan_m3b_m3a[0]
        tan_m3b_m3a_t[1] = tan_m3b_m3a[1] * np.cos(th1_l3) - tan_m3b_m3a[2] * np.sin(
            th1_l3
        )
        tan_m3b_m3a_t[2] = tan_m3b_m3a[1] * np.sin(th1_l3) + tan_m3b_m3a[2] * np.cos(
            th1_l3
        )

        tan_og_t[0] = tan_og_vac[0]
        tan_og_t[1] = tan_og_vac[1] * np.cos(th1_l3) - tan_og_vac[2] * np.sin(th1_l3)
        tan_og_t[2] = tan_og_vac[1] * np.sin(th1_l3) + tan_og_vac[2] * np.cos(th1_l3)

        alpha = tan_og_t[0]
        beta = tan_og_t[1]
        gamma = tan_og_t[2]

        def root_z2b(t):
            x = P_m3a[0] + alpha * t
            y = P_m3a[1] + beta * t
            z = P_m3a[2] + gamma * t

            xm2b, ym2b, zm2b = tele_into_m2b(
                x, y, z
            )  # take ray end coordinates and convert to M1 coordinates
            z_m2b = z2b(xm2b, ym2b)  # Z of mirror 1 in M1 coordinates
            if np.isnan(z_m2b):
                z_m2b = 0
            root = zm2b - z_m2b
            return root

        t_m2b = optimize.brentq(root_z2b, 10, 600)

        # Location of where ray hits M1
        x_m2b = P_m3a[0] + alpha * t_m2b
        y_m2b = P_m3a[1] + beta * t_m2b
        z_m2b = P_m3a[2] + gamma * t_m2b
        P_m2b = np.array([x_m2b, y_m2b, z_m2b])

        ###### in M1 cordinates ##########################
        x_m2b_temp, y_m2b_temp, z_m2b_temp = tele_into_m2b(
            x_m2b, y_m2b, z_m2b
        )  # P_m2a temp
        x_m3a_temp, y_m3a_temp, z_m3a_temp = tele_into_m2b(
            P_m3a[0], P_m3a[1], P_m3a[2]
        )  # P_1b temp
        norm = d_z2b(x_m2b_temp, y_m2b_temp)
        norm_temp = np.array([-norm[0], -norm[1], 1])
        N_hat = norm_temp / np.sqrt(sum(norm_temp ** 2))
        vec_m3a_m2b = np.array([x_m2b_temp, y_m2b_temp, z_m2b_temp]) - np.array(
            [x_m3a_temp, y_m3a_temp, z_m3a_temp]
        )
        dist_m3a_m2b = np.sqrt(np.sum(vec_m3a_m2b ** 2))
        tan_m3a_m2b = vec_m3a_m2b / dist_m3a_m2b

        # Outgoing ray
        tan_og_si = snell_vec(n_vac, n_si, N_hat, tan_m3a_m2b)

        # Transform back to telescope cordinates
        N_hat_t = np.zeros(3)
        tan_m3a_m2b_t = np.zeros(3)
        tan_og_t = np.zeros(3)

        N_hat_t[0] = N_hat[0]
        N_hat_t[1] = N_hat[1] * np.cos(th1_l2) - N_hat[2] * np.sin(th1_l2)
        N_hat_t[2] = N_hat[1] * np.sin(th1_l2) + N_hat[2] * np.cos(th1_l2)

        tan_m3a_m2b_t[0] = tan_m3a_m2b[0]
        tan_m3a_m2b_t[1] = tan_m3a_m2b[1] * np.cos(th1_l2) - tan_m3a_m2b[2] * np.sin(
            th1_l2
        )
        tan_m3a_m2b_t[2] = tan_m3a_m2b[1] * np.sin(th1_l2) + tan_m3a_m2b[2] * np.cos(
            th1_l2
        )

        tan_og_t[0] = tan_og_si[0]
        tan_og_t[1] = tan_og_si[1] * np.cos(th1_l2) - tan_og_si[2] * np.sin(th1_l2)
        tan_og_t[2] = tan_og_si[1] * np.sin(th1_l2) + tan_og_si[2] * np.cos(th1_l2)

        alpha = tan_og_t[0]
        beta = tan_og_t[1]
        gamma = tan_og_t[2]

        def root_z2a(t):
            x = P_m2b[0] + alpha * t
            y = P_m2b[1] + beta * t
            z = P_m2b[2] + gamma * t
            xm2a, ym2a, zm2a = tele_into_m2a(
                x, y, z
            )  # take ray end coordinates and convert to M1 coordinates
            z_m2a = z2a(xm2a, ym2a)  # Z of mirror 1 in M1 coordinates
            if np.isnan(z_m2a):
                z_m2a = 0
            root = zm2a - z_m2a
            return root

        t_m2a = optimize.brentq(root_z2a, 5, 400)

        # Location of where ray hits M1
        x_m2a = P_m2b[0] + alpha * t_m2a
        y_m2a = P_m2b[1] + beta * t_m2a
        z_m2a = P_m2b[2] + gamma * t_m2a
        P_m2a = np.array([x_m2a, y_m2a, z_m2a])

        ###### in M1 cordinates ##########################
        x_m2b_temp, y_m2b_temp, z_m2b_temp = tele_into_m2a(
            x_m2b, y_m2b, z_m2b
        )  # P_m1b temp
        x_m2a_temp, y_m2a_temp, z_m2a_temp = tele_into_m2a(
            P_m2a[0], P_m2a[1], P_m2a[2]
        )  # P_1a temp
        norm = d_z2a(x_m2a_temp, y_m2a_temp)
        norm_temp = np.array([-norm[0], -norm[1], 1])
        N_hat = norm_temp / np.sqrt(sum(norm_temp ** 2))
        vec_m2b_m2a = np.array([x_m2a_temp, y_m2a_temp, z_m2a_temp]) - np.array(
            [x_m2b_temp, y_m2b_temp, z_m2b_temp]
        )
        dist_m2b_m2a = np.sqrt(np.sum(vec_m2b_m2a ** 2))
        tan_m2b_m2a = vec_m2b_m2a / dist_m2b_m2a

        # Outgoing ray

        tan_og_vac = snell_vec(n_si, n_vac, N_hat, tan_m2b_m2a)

        # Transform back to telescope cordinates
        N_hat_t = np.zeros(3)
        tan_m2b_m2a_t = np.zeros(3)
        tan_og_t = np.zeros(3)

        N_hat_t[0] = N_hat[0]
        N_hat_t[1] = N_hat[1] * np.cos(th2_l2) - N_hat[2] * np.sin(th2_l2)
        N_hat_t[2] = N_hat[1] * np.sin(th2_l2) + N_hat[2] * np.cos(th2_l2)

        tan_m2b_m2a_t[0] = tan_m2b_m2a[0]
        tan_m2b_m2a_t[1] = tan_m2b_m2a[1] * np.cos(th2_l2) - tan_m2b_m2a[2] * np.sin(
            th2_l2
        )
        tan_m2b_m2a_t[2] = tan_m2b_m2a[1] * np.sin(th2_l2) + tan_m2b_m2a[2] * np.cos(
            th2_l2
        )

        tan_og_t[0] = tan_og_vac[0]
        tan_og_t[1] = tan_og_vac[1] * np.cos(th2_l2) - tan_og_vac[2] * np.sin(th2_l2)
        tan_og_t[2] = tan_og_vac[1] * np.sin(th2_l2) + tan_og_vac[2] * np.cos(th2_l2)

        alpha = tan_og_t[0]
        beta = tan_og_t[1]
        gamma = tan_og_t[2]

        def root_z1b(t):

            x = x_m2a + alpha * t
            y = y_m2a + beta * t
            z = z_m2a + gamma * t
            xm1b, ym1b, zm1b = tele_into_m1b(
                x, y, z
            )  # Convert ray's endpoint into M2 coordinates
            z_m1b = z1b(xm1b, ym1b)  # Z of mirror in M2 coordinates
            if np.isnan(z_m1b):
                z_m1b = 0
            root = zm1b - z_m1b
            return root

        t_m1b = optimize.brentq(root_z1b, 5, 1000)

        # Location of where ray hits M2
        x_m1b = x_m2a + alpha * t_m1b
        y_m1b = y_m2a + beta * t_m1b
        z_m1b = z_m2a + gamma * t_m1b
        P_m1b = np.array([x_m1b, y_m1b, z_m1b])

        ###### in M2 coordinates ##########################
        x_m1b_temp, y_m1b_temp, z_m1b_temp = tele_into_m1b(
            x_m1b, y_m1b, z_m1b
        )  # P_m2 temp
        x_m2a_temp, y_m2a_temp, z_m2a_temp = tele_into_m1b(
            x_m2a, y_m2a, z_m2a
        )  # P_rx temp
        norm = d_z1b(x_m1b_temp, y_m1b_temp)
        norm_temp = np.array([-norm[0], -norm[1], 1])
        N_hat = norm_temp / np.sqrt(sum(norm_temp ** 2))
        vec_m2a_m1b = np.array([x_m1b_temp, y_m1b_temp, z_m1b_temp]) - np.array(
            [x_m2a_temp, y_m2a_temp, z_m2a_temp]
        )
        dist_m2a_m1b = np.sqrt(np.sum(vec_m2a_m1b ** 2))
        tan_m2a_m1b = vec_m2a_m1b / dist_m2a_m1b

        # Outgoing ray
        tan_og_si = snell_vec(n_vac, n_si, N_hat, tan_m2a_m1b)

        # Outgoing
        # Transform back to telescope cordinates
        N_hat_t = np.zeros(3)
        tan_m2a_m1b_t = np.zeros(3)
        tan_og_t = np.zeros(3)

        N_hat_t[0] = N_hat[0]
        N_hat_t[1] = N_hat[1] * np.cos(th2_l1) - N_hat[2] * np.sin(th2_l1)
        N_hat_t[2] = N_hat[1] * np.sin(th2_l1) + N_hat[2] * np.cos(th2_l1)

        tan_m2a_m1b_t[0] = tan_m2a_m1b[0]
        tan_m2a_m1b_t[1] = tan_m2a_m1b[1] * np.cos(th2_l1) - tan_m2a_m1b[2] * np.sin(
            th2_l1
        )
        tan_m2a_m1b_t[2] = tan_m2a_m1b[1] * np.sin(th2_l1) + tan_m2a_m1b[2] * np.cos(
            th2_l1
        )

        tan_og_t[0] = tan_og_si[0]
        tan_og_t[1] = tan_og_si[1] * np.cos(th2_l1) - tan_og_si[2] * np.sin(th2_l1)
        tan_og_t[2] = tan_og_si[1] * np.sin(th2_l1) + tan_og_si[2] * np.cos(th2_l1)
        ##################################################

        alpha = tan_og_t[0]
        beta = tan_og_t[1]
        gamma = tan_og_t[2]

        def root_z1a(t):
            x = P_m1b[0] + alpha * t
            y = P_m1b[1] + beta * t
            z = P_m1b[2] + gamma * t
            xm1a, ym1a, zm1a = tele_into_m1a(
                x, y, z
            )  # take ray end coordinates and convert to M1 coordinates
            z_m1a = z1a(xm1a, ym1a)  # Z of mirror 1 in M1 coordinates
            if np.isnan(z_m1a):
                z_m1a = 0
            root = zm1a - z_m1a
            return root

        t_m1a = optimize.brentq(root_z1a, 5, 220)

        # Location of where ray hits M1
        x_m1a = P_m1b[0] + alpha * t_m1a
        y_m1a = P_m1b[1] + beta * t_m1a
        z_m1a = P_m1b[2] + gamma * t_m1a
        P_m1a = np.array([x_m1a, y_m1a, z_m1a])

        ###### in M1 cordinates ##########################
        x_m1b_temp, y_m1b_temp, z_m1b_temp = tele_into_m1a(
            x_m1b, y_m1b, z_m1b
        )  # P_m1b temp
        x_m1a_temp, y_m1a_temp, z_m1a_temp = tele_into_m1a(
            P_m1a[0], P_m1a[1], P_m1a[2]
        )  # P_1a temp
        norm = d_z1a(x_m1a_temp, y_m1a_temp)
        norm_temp = np.array([-norm[0], -norm[1], 1])
        N_hat = norm_temp / np.sqrt(sum(norm_temp ** 2))
        vec_m1b_m1a = np.array([x_m1a_temp, y_m1a_temp, z_m1a_temp]) - np.array(
            [x_m1b_temp, y_m1b_temp, z_m1b_temp]
        )
        dist_m1b_m1a = np.sqrt(np.dot(vec_m1b_m1a, vec_m1b_m1a))
        tan_m1b_m1a = vec_m1b_m1a / dist_m1b_m1a

        # Outgoing ray
        tan_og_vac = snell_vec(n_si, n_vac, N_hat, tan_m1b_m1a)

        # Transform back to telescope cordinates
        N_hat_t = np.zeros(3)
        tan_m1b_m1a_t = np.zeros(3)
        tan_og_t = np.zeros(3)

        N_hat_t[0] = N_hat[0]
        N_hat_t[1] = N_hat[1] * np.cos(th1_l1) - N_hat[2] * np.sin(th1_l1)
        N_hat_t[2] = N_hat[1] * np.sin(th1_l1) + N_hat[2] * np.cos(th1_l1)

        tan_m1b_m1a_t[0] = tan_m1b_m1a[0]
        tan_m1b_m1a_t[1] = tan_m1b_m1a[1] * np.cos(th1_l1) - tan_m1b_m1a[2] * np.sin(
            th1_l1
        )
        tan_m1b_m1a_t[2] = tan_m1b_m1a[1] * np.sin(th1_l1) + tan_m1b_m1a[2] * np.cos(
            th1_l1
        )

        tan_og_t[0] = tan_og_vac[0]
        tan_og_t[1] = tan_og_vac[1] * np.cos(th1_l1) - tan_og_vac[2] * np.sin(th1_l1)
        tan_og_t[2] = tan_og_vac[1] * np.sin(th1_l1) + tan_og_vac[2] * np.cos(th1_l1)

        #         ################################################

        dist_m1a_lyot = abs((y_lyot - P_m1a[1]) / tan_og_t[1])
        pos_lyot = P_m1a + dist_m1a_lyot * tan_og_t

        dist_lyot_ap = abs((tele_geo.y_source - pos_lyot[1]) / tan_og_t[1])

        total_path_length = (
            dist_rx_m3b
            + dist_m3b_m3a * (n_si)
            + dist_m3a_m2b
            + dist_m2b_m2a * (n_si)
            + dist_m2a_m1b
            + dist_m1b_m1a * (n_si)
            + dist_m1a_lyot
            + dist_lyot_ap
        )

        pos_ap = pos_lyot + dist_lyot_ap * tan_og_t

        # Estimate theta
        de_ve = np.arctan(tan_rx_m3b_t[0] / (-tan_rx_m3b_t[1]))
        de_ho = np.arctan(
            tan_rx_m3b_t[2] / np.sqrt(tan_rx_m3b_t[0] ** 2 + tan_rx_m3b_t[1] ** 2)
        )

        ################################################

        if ((pos_lyot[0] ** 2 + pos_lyot[2] ** 2) <= 210 ** 2) and (
            (x_m2a ** 2 + z_m2a ** 2) <= (210 ** 2)
        ):

            if plot:
                if np.mod(ii, 53) == 0:
                    alph = 0.2
                    plt.plot([y_0, y_m3b], [z_0, z_m3b], "-", color=col, alpha=alph)
                    plt.plot([y_m3b, y_m3a], [z_m3b, z_m3a], "-", color=col, alpha=alph)
                    plt.plot([y_m3a, y_m2b], [z_m3a, z_m2b], "-", color=col, alpha=alph)
                    plt.plot([y_m2b, y_m2a], [z_m2b, z_m2a], "-", color=col, alpha=alph)
                    plt.plot([y_m2a, y_m1b], [z_m2a, z_m1b], "-", color=col, alpha=alph)
                    plt.plot([y_m1b, y_m1a], [z_m1b, z_m1a], "-", color=col, alpha=alph)
                    plt.plot(
                        [y_m1a, pos_lyot[1]],
                        [z_m1a, pos_lyot[2]],
                        "-",
                        color=col,
                        alpha=alph,
                    )
                    plt.plot(
                        [pos_lyot[1], pos_ap[1]],
                        [pos_lyot[2], pos_ap[2]],
                        "-",
                        color=col,
                        alpha=alph,
                    )

        #         Write out
        out[0, ii] = pos_ap[0]
        out[1, ii] = pos_ap[1]
        out[2, ii] = pos_ap[2]

        if ((pos_lyot[0] ** 2 + pos_lyot[2] ** 2) <= (420 / 2) ** 2) and (
            (x_m2a ** 2 + z_m2a ** 2) <= (210 ** 2)
        ):
            out[3, ii] = total_path_length
            out[4, ii] = np.exp(
                (-0.5)
                * ((de_ho / horn_fwhp_x) ** 2 + (de_ve / horn_fwhp_y) ** 2)
                / (1 / (np.sqrt(8 * np.log(2)))) ** 2
            )
        else:
            out[3, ii] = 0
            out[4, ii] = 0

        out[5, ii] = N_hat_t[0]
        out[6, ii] = N_hat_t[1]
        out[7, ii] = N_hat_t[2]

        out[8, ii] = tan_og_t[0]
        out[9, ii] = tan_og_t[1]
        out[10, ii] = tan_og_t[2]
        pbar.update(1)
    pbar.close()
    return out


def getNearField(tele_geo, rx, plot=False):
    """Get the near field of a receiver feed."""
    # get the ray trace
    out = rx_to_lyot(rx, tele_geo, plot, "b")
    # get the near field
    xx = np.where(out[4] != 0)
    a_sim = out[4][xx]
    p_sim = np.mod(tele_geo.k * (out[3][xx] - np.mean(out[3][xx])) / 1e3 / 2, 2 * np.pi)
    x_sim = out[0][xx] / 1e1
    y_sim = out[2][xx] / 1e1
    tan_og = np.array([np.mean(out[8][xx]), np.mean(out[9][xx]), np.mean(out[10][xx])])

    class SimOutput:
        def __init__(self, x_sim, y_sim, a_sim, p_sim, tan_og):
            self.a_sim = a_sim
            self.p_sim = p_sim
            self.x_sim = x_sim
            self.y_sim = y_sim
            self.tan_og = tan_og

    return SimOutput(x_sim, y_sim, a_sim, p_sim - np.mean(p_sim), tan_og)


def plotSimFields(sb, tele_geo):
    fig, ax = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    ax[1].set_title("Phase")
    cols = ax[1].scatter(sb.x_sim, sb.y_sim, s=5, c=sb.p_sim)
    ax[1].set_aspect("equal")
    ax[1].set_xlim(-40, 40)
    ax[1].set_ylim(-40, 40)
    ax[1].set_xlabel("x (cm)")
    plt.colorbar(cols, ax=ax[1], label="[rad]")
    ax[1].grid()

    ax[0].set_title("Power")
    ax[0].plot(0, 0, "o", color=red)
    cols = ax[0].scatter(
        sb.x_sim,
        sb.y_sim,
        s=2,
        c=20 * np.log10(sb.a_sim / np.max(sb.a_sim)),
        vmax=0,
        vmin=-10,
    )
    ax[0].set_aspect("equal")
    ax[0].set_xlim(-40, 40)
    ax[0].set_ylim(-40, 40)
    ax[0].set_xlabel("x (cm)")
    ax[0].set_ylabel("y (cm)")
    ax[0].grid()
    plt.colorbar(cols, ax=ax[0], label="dB")
    plt.show()


def get_sim(sb):
    return sb.a_sim, sb.p_sim, sb.x_sim, sb.y_sim


def get_beam_cent(sb):
    beam_cent = [np.mean(sb.x_sim), np.mean(sb.y_sim)]
    return beam_cent


def get_fwhm(sb):
    a_sim, p_sim, x_sim, y_sim = get_sim(sb)
    beam_cent = get_beam_cent(sb)
    fwhm_nf_x = abs(y_sim - beam_cent[1])[
        np.where(a_sim ** 2 > np.max(a_sim ** 2) / 2)
    ].max()
    fwhm_nf_y = abs(x_sim - beam_cent[0])[
        np.where(a_sim ** 2 > np.max(a_sim ** 2) / 2)
    ].max()
    return fwhm_nf_x, fwhm_nf_y


def get_angle_out(sb):
    theta_x = abs(90 - np.rad2deg(np.arctan2(sb.tan_og[1], sb.tan_og[2])))
    theta_y = abs(90 - np.rad2deg(np.arctan2(sb.tan_og[1], sb.tan_og[0])))
    return [theta_x, theta_y]


def sim2d(sb):

    x_sim = sb.x_sim  # [cm]
    y_sim = sb.y_sim  # [cm]

    beam_cent = get_beam_cent(sb)
    indx_xy = np.where(
        (np.isnan(x_sim) == False)
        & (np.isnan(y_sim) == False)
        & ((x_sim - beam_cent[0]) ** 2 + (y_sim - beam_cent[1]) ** 2 <= 21.0 ** 2)
    )
    grid_x, grid_y = np.mgrid[
        np.min(x_sim[indx_xy]) : np.max(x_sim[indx_xy]) : 100j,
        np.min(y_sim[indx_xy]) : np.max(y_sim[indx_xy]) : 100j,
    ]
    points = np.array([x_sim[indx_xy], y_sim[indx_xy]]).T
    values = np.array(sb.p_sim[indx_xy]).T
    from scipy.interpolate import griddata

    grid_z2 = griddata(points, values, (grid_x, grid_y), method="cubic")
    grid_z2 = np.where(
        ((grid_x - (beam_cent[0])) ** 2 + (grid_y - (beam_cent[1])) ** 2 >= 20 ** 2),
        0,
        grid_z2,
    )
    grid_z2 = np.where(np.isnan(grid_z2), 0, grid_z2)
    values = np.array(sb.a_sim[indx_xy]).T
    grid_amp = griddata(points, values, (grid_x, grid_y), method="cubic")
    grid_amp = np.where(
        ((grid_x - (beam_cent[0])) ** 2 + (grid_y - (beam_cent[1])) ** 2 >= 20 ** 2),
        0,
        grid_amp,
    )
    grid_amp = np.where(np.isnan(grid_amp), 0, grid_amp)

    class SimOutput:
        def __init__(self, x_sim, y_sim, a_sim, p_sim):
            self.a_sim = grid_amp
            self.p_sim = grid_z2
            self.x_sim = grid_x
            self.y_sim = grid_y

    return SimOutput(grid_x, grid_y, grid_amp, grid_z2)


def sim_data(sb_in, STAGE_RANGE, STEP):
    """
    STAGE_RANGE : range of stage movement in cm
    STEP : step size in cm
    """
    beam_final = sb_in.a_sim.T * np.exp(complex(0, 1) * sb_in.p_sim.T)

    x_interp = sb_in.x_sim[:, 0]
    y_interp = sb_in.y_sim[0, :]
    func_beam = interpolate.interp2d(
        x_interp,
        y_interp,
        np.abs(beam_final) / np.max(np.abs(beam_final)),
        kind="cubic",
    )
    func_phase = interpolate.interp2d(
        x_interp,
        y_interp,
        np.arctan2(np.imag(beam_final), np.real(beam_final)),
        kind="cubic",
    )
    x_data = np.arange(-STAGE_RANGE / 2, STAGE_RANGE / 2, STEP)
    y_data = np.arange(-STAGE_RANGE / 2, STAGE_RANGE / 2, STEP)
    beam_data = func_beam(x_data, y_data)
    phase_data = func_phase(x_data, y_data)
    x_data, y_data = np.meshgrid(x_data, y_data)

    class output:
        def __init__(self, x_sim, y_sim, a_sim, p_sim):
            self.a_data = beam_data / np.max(beam_data)
            self.p_data = phase_data
            self.x_data = x_data
            self.y_data = y_data

    return output(x_data, y_data, beam_data, phase_data)


def ff_fwhm_est(beam_fft, x_ang, y_ang):

    a = abs(beam_fft) ** 2 / np.max(abs(beam_fft) ** 2)
    x_out = y_ang
    y_out = x_ang

    indx = np.where(abs(a) == np.max(abs(a)))
    x = y_out[indx[0][0], :] * 60 * 180 / np.pi
    y = abs(a)[indx[0][0], :] / np.max(abs(a))
    v1 = x[np.where((y > 0.5))][0]
    v2 = x[np.where((y > 0.5))][-1]

    fwhm1 = abs(v1 - v2)

    indx = np.where(abs(a) == np.max(abs(a)))
    x = x_out[:, indx[1][0]] * 60 * 180 / np.pi
    y = abs(a)[:, indx[1][0]] / np.max(abs(a))
    v1 = x[np.where((y > 0.5))][0]
    v2 = x[np.where((y > 0.5))][-1]

    fwhm2 = abs(v1 - v2)
    fwhm = (fwhm1 + fwhm2) / 2
    return fwhm
