import numpy as np
import os
from ...helpers.misc import import_config


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
