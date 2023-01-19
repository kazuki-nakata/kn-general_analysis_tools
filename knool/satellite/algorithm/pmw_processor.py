import numpy as np
import os
from ...helpers.misc import import_config
from scipy.special import jv
from ...geodata import geo_info
from ...fortlib import pmw_processor as fpmw


def calc_TBD_TBU_Tau(wv, ts, clw, eaz, freq_list=["6.9GHz", "18.7GHz", "23.8GHz", "36.5GHz", "89.0GHz"]):
    # bulk formula from Wentz (2000) simple RTM
    config_path = os.path.dirname(__file__) + os.sep + "conf/amsr_wentz.yaml"
    params = import_config(config_path=config_path)

    rad = np.radians(eaz)
    tv = np.where(wv <= 48, 273.16 + 0.8337 * wv - 3.029 * (10 ** (-5)) * wv**3.33, 301.16)
    zeta = np.where(np.abs(ts - tv) <= 20, 1.05 * (ts - tv) * (1 - (ts - tv) ** 2 / 1200), np.sign(ts - tv) * 14)

    output = []

    for freq in freq_list:
        p = params["atmosphere"][freq]
        td = p["b0"] + p["b1"] * wv + p["b2"] * wv**2 + p["b3"] * wv**3 + p["b4"] * wv**4 + p["b5"] * zeta
        tu = td + p["b6"] + p["b7"] * wv

        ao = p["aO1"] + p["aO2"] * (td - 270)
        av = p["aV1"] * wv + p["aV2"] * wv**2
        al = p["aL1"] * (1 - p["aL2"] * ((ts + 273) / 2 - 283)) * clw
        tau = np.exp(-(ao + av + al) / np.cos(rad))
        alpha = 1 - tau
        result = np.array([alpha * td, alpha * tu, tau])
        output.append(result)

    return np.array(output)


def calc_emissivity_from_SRTM(tb, tbd, tbu, tau, sst, tbc=2.7):
    # Calclate mean emissivity from measured (Tb, sst) and modelled (tbd,tbu,tau) from calc_TBD_TBU_Tau function.
    em = (tb - tbu - tau * tbd - tau * tau * tbc) / tau / (sst - tbd)
    return em


def calc_tb_from_SRTM(em, sst, tbd, tbu, tau, tbc=2.7):
    # calculate mean brightness temperature from emissivity, sst, tbd, tbu, and tau
    return tbu + tau * (em * sst + (1 - em) * tbd + tau * tbc)


def get_antenna_pattern_gaussian_beam(fwhm_x, fwhm_y, x_array, y_array, offset_x=0, offset_y=0):
    # fwhm:footprint size(-3db) nx,ny: km for azimuth and elevation axis
    sigma_x = fwhm_x / np.sqrt(2 * np.log(2))
    sigma_y = fwhm_y / np.sqrt(2 * np.log(2))
    Z = np.exp(-2 * ((x_array - offset_x) ** 2 / (sigma_x**2) + (y_array - offset_y) ** 2 / (sigma_y**2)))
    return Z


def get_antenna_pattern_bessel_beam(fwhm_x, fwhm_y, x_array, y_array, offset_x=0, offset_y=0):
    # fwhm:footprint size(-3db) nx,ny: km for azimuth and elevation axis
    X2 = (x_array - offset_x) ** 2
    Y2 = (y_array - offset_y) ** 2
    s = 0.5 / 3.2106
    r = 1 / s * np.sqrt(X2 / (fwhm_x**2) + Y2 / (fwhm_y**2))
    j = jv(3, r)
    Z = 47.9985 * j / (r**3)
    return Z


def antenna_pattern_integration(integ_time, integ_interval, antenna_func, func_args, rot_velo=40.0):

    itime_list = np.arange(-integ_time / 2, integ_time / 2 + integ_interval, integ_interval) * 10 ** (-3)  # second
    r = rot_velo * 2 * np.pi / 60  # radian/s
    ap = 0
    for itime in itime_list:
        offset_x = r * itime * 180 / 3.14
        az = np.arange(-3, 3, 0.01)
        el = np.arange(-3, 3, 0.01)
        X, Y = np.meshgrid(az, el)
        ap0 = antenna_func(*func_args, offset_x=offset_x)
        ap = ap + ap0
    ap = ap / np.max(ap)
    return ap


def run_rSIR(grid_x, grid_y, eaz, time, tb, mask, lon, lat, day, ltod_min, ltod_max, rot_offset, wsize, res):
    test, test2 = fpmw.rsir(
        grid_x,
        grid_y,
        eaz,
        time,
        tb,
        mask,
        lon,
        lat,
        int(day),
        float(ltod_min),
        float(ltod_max),
        float(rot_offset),
        int(wsize),
        float(res),
    )
