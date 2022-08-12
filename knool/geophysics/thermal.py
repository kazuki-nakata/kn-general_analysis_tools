import os
import numpy as np
from ..helpers.misc import import_config

print(os.path.dirname(__file__) + os.sep + "params_physics_aoi.yaml")
params = import_config(config_path=os.path.dirname(__file__) + os.sep + "params_physics_aoi.yaml")


def calc_ice_production(hb):
    dt = 24.0 * 60.0 * 60.0
    dhi = np.abs(hb) * dt / params["seaice"]["rho"] / params["seaice"]["lf"]
    return dhi


# def convert_td_to_q_K2002(td2m): #Kspecific humidity ara et al. (2002)


def calc_vapor_pressure(t, mode=1):  # input temperature (Kelvin)
    if mode == 1:  # Kondo 1975 (described in MRI manual) for ocean
        t = t - 273.15
        ea = 0.98 * 6.1078 * 10 ** (7.5 * t / (273.3 + t))
    elif mode == 2:  # Ohshima et al. 2003 for ice and ocean
        t = t - 273.15
        ea = 6.11 * 10 ** (9.5 * t / (265.3 + t))
    elif mode == 3:  # Yin for ocean
        t = t - 273.15
        ea = 6.11 * 10 ** (7.5 * t / (-35.86 + t + 273.15))
        # ea = 6.1078 * np.exp(17.269 * t / (237.3 + t)) #almost the same as the above fomula
    elif mode == 4:  # LV1990 for any surface?
        ea = np.where(t > 273.15, np.exp((-6763.6 / t) - 4.9283 * np.log(t) + 54.23), np.exp((-6141.0 / t) + 24.3))
    return ea


def calc_swave(al, clo, td2m, jday, lat, hour):
    s = 1358.0
    b01 = 9.5
    b02 = 265.3
    a1 = -23.44
    a2 = 172.0
    st = 0.0
    ha = 0.0
    sz = 0.0
    rjday = jday

    #   --- calc the vapor pressure -----------------------------
    ev = 6.11 * 10.0 ** ((b01 * (td2m - 273.15)) / (b02 + td2m - 273.15))
    #   --- calc dec: declination -------------------------------
    dec = np.radians(a1 * np.cos(np.radians(a2 - rjday)))
    rad_la = np.radians(lat)

    st = hour
    ha = (12.0 - st) * np.pi / 12.0
    sz = np.sin(rad_la) * np.sin(dec) + np.cos(rad_la) * np.cos(dec) * np.cos(ha)
    if sz <= 0.0:
        sz = 0.0
    q1 = np.where(sz < 0, 0, (s * (sz**2)) / ((sz + 2.7) * ev * (1.0e-3) + 1.085 * sz + 0.10))

    an = np.angle(np.arcsin(sz))
    q2 = (1 - al) * q1 * (1 - 0.62 * clo + 0.0019 * an)  # (Andreas and Ackley,1982)
    return q2


def calc_lwave_MC1973(clo, st, t2m, em):  # maykut and churtch 1973 by alaska data
    # --- incoming longwave radiation ---
    ila = (
        0.7855 * (1.0 + 0.2232 * clo**2.75) * params["general"]["sb_const"] * t2m**4.0
    )  # (Maykut and Church, 1973)
    ola = -(em * params["general"]["sb_const"] * st**4)
    return ila, ola


def calc_lwave_KA1994(clo, st, t2m, em):  # Koenig-Langlo and Augstein, 1994 by polar region data
    ila = (0.765 + 0.22 * clo**3.0) * params["general"]["sb_const"] * t2m**4
    ola = -(em * params["general"]["sb_const"] * st**4)
    return ila, ola


def calc_lwave_G1998(clo, st, t2m, em):  # Guest（1997）by weddell data
    ila = (params["general"]["sb_const"] * t2m**4 - 85.6) * (1 + 0.26 * clo)
    ola = -(em * params["general"]["sb_const"] * st**4)
    return ila, ola


def calc_lwave_J2006(st, t2m, td2m, em):  # Jun et al. (2006) for arctic
    # ea = (qa * p) / 0.622
    ea = calc_vapor_pressure(td2m, mode=4)
    eps_a = (0.0003 * (t2m - 273.16) ** 2 - 0.0079 * (t2m - 273.16) + 1.2983) * (ea / t2m) ** (1 / 7)
    ila = eps_a * params["general"]["sb_const"] * t2m**4
    ola = -(em * params["general"]["sb_const"] * st**4)
    return ila, ola


def calc_transfer_coeffs_K1975(wg, t2m, tsfc):

    params2 = import_config(config_path=os.path.dirname(__file__) + os.sep + "params_heat.yaml")
    params2 = params2["Kondo1975"]

    # --- set PARAMETER ---
    # ---- bulk transfer coefficients ---
    calc_bulk_cha = lambda wg, p: (p["ah"] + (p["bh"] * (wg ** p["ph"])) + (p["ch"] * ((wg - 8) ** 2))) / 1000.0
    calc_bulk_cea = lambda wg, p: (p["ae"] + (p["be"] * (wg ** p["pe"])) + (p["ce"] * ((wg - 8) ** 2))) / 1000.0

    cha = np.where(
        wg < 2.2,
        calc_bulk_cha(wg, params2[0]),
        np.where(
            (wg >= 0.02) & (wg < 5.0),
            calc_bulk_cha(wg, params2[1]),
            np.where(
                (wg >= 5.0) & (wg < 8.0),
                calc_bulk_cha(wg, params2[2]),
                np.where(
                    (wg >= 8.0) & (wg < 25.0),
                    calc_bulk_cha(wg, params2[3]),
                    np.where((wg > 25.0) & (wg < 50.0), calc_bulk_cha(wg, params2[4]), np.NaN),
                ),
            ),
        ),
    )

    cea = np.where(
        wg < 2.2,
        calc_bulk_cea(wg, params2[0]),
        np.where(
            (wg >= 0.02) & (wg < 5.0),
            calc_bulk_cea(wg, params2[1]),
            np.where(
                (wg >= 5.0) & (wg < 8.0),
                calc_bulk_cea(wg, params2[2]),
                np.where(
                    (wg >= 8.0) & (wg < 25.0),
                    calc_bulk_cea(wg, params2[3]),
                    np.where((wg > 25.0) & (wg < 50.0), calc_bulk_cea(wg, params2[4]), np.NaN),
                ),
            ),
        ),
    )

    # --- 安定度 ---
    dtsfc = tsfc - t2m
    s0 = dtsfc * wg ** (-2)
    s = s0 * (np.abs(s0) / (np.abs(s0) + 0.01))

    # ---calc transfer coefficients----------------------
    ch = np.where(
        dtsfc == 0,
        cha,  # 中立
        np.where(
            dtsfc > 0,
            cha * (1.0 + 0.63 * np.sqrt(np.abs(s))),  # 不安定。np.absはwarning出力を抑えるため。
            np.where((s > -3.3) & (s < 0.0), cha * (0.1 + 0.03 * s + 0.9 * np.exp(4.8 * s)), 0),
        ),
    )  # 安定

    ce = np.where(
        dtsfc == 0,
        cea,  # 中立
        np.where(
            dtsfc > 0,
            cea * (1.0 + 0.63 * np.sqrt(np.abs(s))),  # 不安定
            np.where((s > -3.3) & (s < 0.0), cea * (0.1 + 0.03 * s + 0.9 * np.exp(4.8 * s)), 0),
        ),
    )  # 安定

    return ch, ce


def calc_theat_O2003(ti, tw, wg, ic, t2m, td2m, slp):
    wg = np.where(wg == 0, 0.1, wg)
    tsfc = np.where(ic < 0.15, tw, ic * ti + (1.0 - ic) * tw)

    ch, ce = calc_transfer_coeffs_K1975(wg, t2m, tsfc)

    b01 = 9.5
    b02 = 265.3
    rhoa = params["air"]["rho"]
    cp = params["air"]["cp"]
    # ! --- Sensible heat ----------------------------------------
    sehi = rhoa * cp * ch * wg * (t2m - ti)
    sehw = rhoa * cp * ch * wg * (t2m - tw)
    # ! --- Latent heat ------------------------------------------
    ea = 6.11 * 10.0 ** ((b01 * (td2m - 273.15)) / (b02 + td2m - 273.15))
    esw = 6.11 * 10.0 ** ((b01 * (tw - 273.15)) / (b02 + tw - 273.15))
    esi = 6.11 * 10.0 ** ((b01 * (ti - 273.15)) / (b02 + ti - 273.15))
    lahw = 0.622 * rhoa * 2.52e6 * ce / (slp / 100.0) * wg * (ea - esw)
    lahi = (
        0.622 * rhoa * 2.86e6 * ce / (slp / 100.0) * wg * (ea - esi)
    )  # latent heat of fusion of ice plus evapolation, thus latent heat for sublimation
    return sehi, sehw, lahi, lahw


# fmt: off
def convert_thermal_ice_thickness(q, ti, hs, tb):
    cki = params["seaice"]["ck"]
    cks = params["snow"]["ck"]
    hi = np.where(q > 0, 1.0, (cki * cks / q * (ti - tb) - cki * hs) / cks)
    return hi
# fmt: on


# fmt: off
def calc_thermal_ice_properties(slp, t2m, td2m, w10m, ic, cl, tw, lat, jday, hour, hs, ti,
    lw_func=calc_lwave_MC1973,
    sw_func=None,
    th_func=calc_theat_O2003,
    test=False):

    emi = params["seaice"]["em"]
    emw = params["seawater"]["em"]
    ali = params["seaice"]["albedo_10_20cm"]
    alw = params["seawater"]["albedo"]

    ilw_ice = olw_ice = 0
    ilw_water = olw_water = 0
    sw_ice = sw_water = 0
    sh_ice = sh_water = lh_ice = lh_water = 0

    if lw_func is not None:
        if lw_func.__name__ == "calc_lwave_J2006":
            ilw_ice, olw_ice = lw_func(ti, t2m, td2m, emi)
            ilw_water, olw_water = lw_func(ti, t2m, td2m, emw)
        else:
            ilw_ice, olw_ice = lw_func(cl, ti, t2m, emi)
            ilw_water, olw_water = lw_func(cl, tw, t2m, emw)
    if sw_func is not None:
        sw_ice = sw_func(ali, cl, td2m, jday, lat, hour)
        sw_water = sw_func(alw, cl, td2m, jday, lat, hour)
    if th_func is not None:
        sh_ice, sh_water, lh_ice, lh_water = th_func(ti, tw, w10m, ic, t2m, td2m, slp)
    hbi_all = (sw_ice + ilw_ice + olw_ice + sh_ice + lh_ice) * ic
    hbw_all = (sw_water + ilw_water + olw_water + sh_water + lh_water) * (1.0 - ic)
    hb_all = hbw_all + hbi_all

    tb = 273.15 - 1.86
    ice_th = convert_thermal_ice_thickness(hb_all, ti, hs, tb)

    if test:
        hbi_list = [sw_ice, ilw_ice, olw_ice, sh_ice, lh_ice]
        hbw_list = [sw_water, ilw_water, olw_water, sh_water, lh_water]
        return hbi_all, hbw_all, ice_th, hbi_list, hbw_list
    else:
        return hbi_all, hbw_all, ice_th
