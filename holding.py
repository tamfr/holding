# coding: utf-8
# Python 3.6
# Holding entries via parametric equations

import numpy as np
import matplotlib.pyplot as plt


def normalize_deg(deg_to_norm):
    """ Ensures degrees are between 0 and 360. """
    add_360 = (deg_to_norm < 0) * 360
    sub_360 = (deg_to_norm > 360) * 360
    return deg_to_norm + add_360 - sub_360


def plot_entry_orbit(
        hdg,
        inbound_crs,
        R,
        v,
        turn_rate,
        wind_dir,
        wind_vel):

    adjusted_hdg = normalize_deg(hdg - inbound_crs + 180)

    switch_sign = (-1)**R

    if 0 + 291*R <= adjusted_hdg <= 69 + 291*R:
        entry_type = 'Teardrop Entry'
        teardrop_hdg = normalize_deg(inbound_crs + 180 + 30 * switch_sign)

        if hdg > teardrop_hdg:
            turn = hdg - (180 + 30*switch_sign) - inbound_crs
            R = (not R)*1*R
        else:
            turn = 180 + 30*switch_sign + inbound_crs - hdg
            R = (not R)*1 + R

    elif 69 + 40*R < adjusted_hdg < 251 + 40*R:
        entry_type = 'Direct Entry'
        turn = (hdg - inbound_crs)*switch_sign + 180

    else:
        entry_type = 'Parallel Entry'
        turn = (180 - hdg)*switch_sign + inbound_crs
        R = (not R)*1

    hdg_xform = - hdg * np.pi / 180
    angular_vel = turn_rate * np.pi / 180  # rad/s
    wind_xform = - wind_dir * np.pi / 180 + 3 * np.pi / 2  # rad

    r = v / 3600 / angular_vel
    turn = normalize_deg(turn)
    time_to_turn = (turn / turn_rate)
    t = np.arange(0, time_to_turn, np.pi / 50)

    wind_x = wind_vel / 3600 * np.cos(wind_xform)
    wind_y = wind_vel / 3600 * np.sin(wind_xform)

    phase = hdg_xform + R * np.pi

    # Shift starting point to origin
    origin_x = r * np.cos(phase)
    origin_y = r * np.sin(phase)

    theta = angular_vel * t * (-1)**R

    x = r * np.cos(theta + phase) + wind_x * t - origin_x
    y = r * np.sin(theta + phase) + wind_y * t - origin_y

    inbound_crs_xform = (180 - inbound_crs) * np.pi / 180 + np.pi / 2

    r = np.arange(0, 4, 0.1)

    x_crs = r*np.cos(inbound_crs_xform)
    y_crs = r*np.sin(inbound_crs_xform)

    ax = plt.subplot(111, polar=False)
    ax.plot(x, y)
    ax.plot(x_crs, y_crs)
    ax.grid(True)
    ax.set(aspect=1, adjustable='datalim')
    ax.set_title(entry_type, va='bottom')

    return plt.show()


def plot_holding_basic_area():

    a_l = 4.5
    l_m = 4.3
    m_g = 5.6
    l_i = 3.5  # Also M-H
    m_e = 5.3
    a_b = 1.5  # Also G-F and J-K

    l_b = (a_l**2 + a_b**2)**(1 / 2)

    step = 0.01

    x_a_l = np.arange(- a_l, 0, step)
    y_a_l = x_a_l * 0
    x_l_m = np.arange(0, l_m, step)
    y_l_m = x_l_m * 0
    x_m_g = np.arange(l_m, m_g + l_m, step)
    y_m_g = x_m_g * 0

    # Calculate arc CB
    l_e = (l_m ** 2 + m_e ** 2) ** (1 / 2)
    phi = np.arccos(l_b / l_e)
    beta = np.arctan(m_e / l_m)
    epsilon = phi + beta
    zeta = np.pi - np.arccos(a_l / l_b)
    theta = np.arange(epsilon, zeta + step, step)
    x_arc_c_b = l_b * np.cos(theta)
    y_arc_c_b = l_b * np.sin(theta)

    # Calculate arc BI
    focal = calc_focal_point(l_i, a_l, a_b)

    epsilon = np.arctan((a_b - focal[1]) / (-a_l - focal[0]))
    zeta = np.arctan((-focal[0]) / (-l_i - focal[1]))

    theta = np.arange(np.pi + epsilon, 3*np.pi/2 - zeta + step, step)
    x_arc_b_i = focal[2] * np.cos(theta) + focal[0]
    y_arc_b_i = focal[2] * np.sin(theta) + focal[1]

    # Calculate arc HG
    focal = calc_focal_point(l_i, m_g, a_b)
    epsilon = np.arctan(focal[0] / (focal[1] + l_i))
    zeta = np.arctan((focal[1] - a_b) / (focal[0] + m_g))
    theta = np.arange(3 * np.pi/2 + epsilon, 2 * np.pi - zeta + step, step)
    x_arc_h_g = focal[2] * np.cos(theta) - focal[0] + l_m
    y_arc_h_g = focal[2] * np.sin(theta) + focal[1]

    # Calculate arc FE
    focal = calc_focal_point(m_e, m_g, -a_b)
    epsilon = np.arctan((focal[1] + a_b) / (m_g - focal[0]))
    zeta = np.arctan(focal[0] / (focal[1] + m_e))
    theta = np.arange(epsilon, np.pi/2 + zeta + step, step)
    x_arc_f_e = focal[2] * np.cos(theta) + focal[0] + l_m
    y_arc_f_e = focal[2] * np.sin(theta) - focal[1]

    x_perimeter = np.concatenate([x_arc_f_e, x_arc_c_b, x_arc_b_i, x_arc_h_g])
    y_perimeter = np.concatenate([y_arc_f_e, y_arc_c_b, y_arc_b_i, y_arc_h_g])

    ax = plt.subplot(111, polar=False)

    ax.plot(x_a_l, y_a_l)
    ax.plot(x_l_m, y_l_m)
    ax.plot(x_m_g, y_m_g)

    ax.plot(x_perimeter, y_perimeter)

    ax.grid(True)
    ax.set(aspect=1, adjustable='datalim')

    return plt.show()


def calc_focal_point(f, lateral_offset, vertical_offset):
    """ Triangulates the focal point of the arc. """
    r = (vertical_offset**2 + lateral_offset**2)**(1/2)
    h = ((f + vertical_offset)**2 + lateral_offset**2)**(1/2)
    omega = np.arccos(h / (2*r))
    phi = np.arccos((h**2 + r**2 - f**2) / (2*h*r))
    gamma = omega - phi
    c = (2*r**2 * (1 - np.cos(gamma)))**(1/2)
    psi = np.arccos((f**2 + c**2 - r**2) / (2*f*c))
    mu = psi - np.pi/2

    x_focal = c * np.cos(mu)
    y_focal = c * np.sin(mu)

    return x_focal, y_focal, r


heading = 350  # degrees
right_turn = 1
inbound_course = 270  # degrees

velocity = 250  # KIAS
rate_of_turn = 3  # degrees/s

wind_velocity = 20  # kts
wind_direction = 180  # degrees

# entry_orbit = plot_entry_orbit(
#     heading,
#     inbound_course,
#     right_turn,
#     velocity,
#     rate_of_turn,
#     wind_direction,
#     wind_velocity)

plot_holding_basic_area()
