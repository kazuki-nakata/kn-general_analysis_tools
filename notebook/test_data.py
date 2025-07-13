import os
import sys
import numpy as np
import scipy
from scipy import signal


def make_wave_3D(frequency, prop_direction, dispersion, size, init_phase=0):
    frequency = frequency*np.pi/180
    prop_direction = prop_direction*np.pi/180
    wave_number = dispersion*frequency  # dispersion_relation
    init_phase = init_phase*np.pi/180
#    phase_velocity=frequency/wave_number
    kx = wave_number*np.cos(prop_direction)
    ky = wave_number*np.sin(prop_direction)
    tmax, ymax, xmax = size

    interval = np.linspace(1, xmax, xmax)
    x0 = np.array([interval for i in range(ymax)])
    x = np.array([x0 for i in range(tmax)])

    interval = np.linspace(1, ymax, ymax)
    y0 = np.array([interval for i in range(xmax)]).T
    y = np.array([y0 for i in range(tmax)])

    interval = np.linspace(1, tmax, tmax)
    t0 = np.array([interval for i in range(ymax)])
    t = np.array([t0 for i in range(xmax)]).T
    output = np.cos(kx*x+ky*y+frequency*t+init_phase) + \
        1j*np.sin(kx*x+ky*y+frequency*t+init_phase)
    return output


def make_stationary_wave_3D(frequency, prop_direction, dispersion, size):
    output1 = make_wave_3D(frequency, prop_direction, dispersion, size)
    output2 = make_wave_3D(frequency, prop_direction+90, dispersion, size)
    return output1+output2
