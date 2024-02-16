import numpy as np
import os
from ...helpers.misc import import_config
from scipy.special import jv
from ...geodata import geo_info
from ...fortlib import pmw_processor as fpmw
from ...fortlib import pmw_processor_l1b as fpmwl1b
from ...fortlib import grid_data
from ...fortlib import sensor_geometry as sg
from ...geodata import geo_map
from osgeo import osr
from ...image import img_destripe


def calc_TBD_TBU_Tau(wv, ts, clw, eaz, freq_list=["6.9GHz", "18.7GHz", "23.8GHz", "36.5GHz", "89.0GHz"]):
    '''
    bulk formula from Wentz (2000) simple RTM
    wv: columnar water vapor (mm)
    clw: cloud liquid water (mm)
    '''
    config_path = os.path.dirname(__file__) + os.sep + "conf/amsr_wentz.yaml"
    params = import_config(config_path=config_path)

    rad = np.radians(eaz)
    tv = np.where(wv <= 48, 273.16 + 0.8337 * wv -
                  3.029 * (10 ** (-5)) * wv**3.33, 301.16)
    zeta = np.where(np.abs(ts - tv) <= 20, 1.05 * (ts - tv) *
                    (1 - (ts - tv) ** 2 / 1200), np.sign(ts - tv) * 14)

    output = []

    for freq in freq_list:
        p = params["atmosphere"][freq]
        td = p["b0"] + p["b1"] * wv + p["b2"] * wv**2 + \
            p["b3"] * wv**3 + p["b4"] * wv**4 + p["b5"] * zeta
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
    # Calclate mean emissivity from measured (Tb, sst) using (tbd,tbu,tau) from calc_TBD_TBU_Tau function.
    em = (tb - tbu - tau * tbd - tau * tau * tbc) / tau / (sst - tbd)
    return em


def calc_tb_from_SRTM(em, sst, tbd, tbu, tau, tbc=2.7):
    # calculate mean brightness temperature from emissivity, sst, tbd, tbu, and tau
    return tbu + tau * (em * sst + (1 - em) * tbd + tau * tbc)


def get_antenna_pattern_gaussian_beam(fwhm_x, fwhm_y, x_array, y_array, offset_x=0, offset_y=0):
    # fwhm:footprint size(-3db) nx,ny: km for azimuth and elevation axis
    sigma_x = fwhm_x / np.sqrt(2 * np.log(2))
    sigma_y = fwhm_y / np.sqrt(2 * np.log(2))
    Z = np.exp(-2 * ((x_array - offset_x) ** 2 / (sigma_x**2) +
               (y_array - offset_y) ** 2 / (sigma_y**2)))
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

    itime_list = np.arange(-integ_time / 2, integ_time / 2 +
                           integ_interval, integ_interval) * 10 ** (-3)  # second
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


def calc_boresight_basis_vectors(p, s):
    # p: ecef obs. location vector at earth surface
    # s: ecef s/c position vector
    i, j, k = sg.calc_boresight_basis_vectors(p, s)
    return i, j, k


def calc_local_az_el_angle(i, j, k, b, p0, p):
    # i,j,k:boresight basis vectors
    # b: boresight vector
    # p0: ecef obs location on earth
    # p: ecef local obs location on earth
    if len(i.shape) == 1:
        az, el = sg.calc_local_az_el_angle(i, j, k, b, p0, p)
    elif len(i.shape) == 2:
        az, el = sg.calc_local_az_el_angle2(i, j, k, b, p0, p)
    return az, el


def l1B_to_l3_WA_TYPE1(grid_x, grid_y, val, mask, wsize, fwhm, res, rm_outer):
    sigma = fwhm / (2 * np.sqrt(np.log(2)))
    output = grid_data.weighted_mean_sigma(
        grid_x, grid_y, val, mask, wsize, sigma, res, rm_outer)
    return output


def l1B_to_l3_WA_TYPE2(grid_x, grid_y, vs, vb, vg, val, mask, wsize, ap, int_ap, res, fwhm, rm_outer):
    sigma = fwhm / (2 * np.sqrt(np.log(2)))
    output = grid_data.weighted_mean_sat(
        grid_x, grid_y, vs, vb, vg, val, mask, wsize, ap, int_ap, res, sigma, rm_outer)
    return output


def run_rSIR(grid_x, grid_y, vs, vb, vg, tbv, mask, wsize, ap, int_ap, res, fwhm, iterate):
    output = fpmw.rsir(grid_x, grid_y, vs, vb, vg, tbv,
                       mask, wsize, ap, int_ap, res, fwhm, iterate)
    return output


def run_rSIRvh(grid_x, grid_y, vs, vb, vg, tb, mask, wsize, ap, int_ap, res, fwhm, iterate):
    output = fpmw.rsirvh(grid_x, grid_y, vs, vb, vg, tb,
                         mask, wsize, ap, int_ap, res, fwhm, iterate)
    return output


def run_rSIR_SS(grid_x, grid_y, grid_id, vs, vb, vg, tbv, mask, id1, id2, wsize, ap, int_ap, res, fwhm, iterate):
    output = fpmw.rsir2(grid_x, grid_y, grid_id, vs, vb, vg, tbv,
                        mask, id1, id2, wsize, ap, int_ap, res, fwhm, iterate)
    return output


def run_re_iteration_l1b(grid_x, grid_y, vs, vb, vg, tbv, mask, init, wsize, ap, int_ap, wr, params, method):
    if method == "rsir":
        output = fpmwl1b.rsir(grid_x, grid_y, vs, vb, vg, tbv,
                              mask, init, wsize, ap, int_ap, wr, params)
    elif method == "banach":
        output = fpmwl1b.banach_gradient(grid_x, grid_y, vs, vb, vg, tbv,
                                         mask, init, wsize, ap, int_ap, wr, params[0], params[1], params[2])
    return output


def interpolation(latb, lonb, imgb, imgb_ij, lata, lona, sp, imga_ij, ap, offset, range_x=5, range_y=5):
    """For How To Use, Prease refer to the explanation of 89GHz interpolation method"""
    # 初期化
    n_proc = 15
    imgb_nx = latb.shape[0]
    imgb_ny = latb.shape[1]
    imgb_a = np.zeros(lata.shape)

    # インデックスカーネルの作成
    cols = np.arange(-range_y, range_y)
    rows = np.arange(-range_x, range_x)
    index_i, index_j = np.meshgrid(cols, rows)
    index_ij = np.array(
        [[index_i.reshape(-1), index_j.reshape(-1)]]).astype(np.int32)
    length = imga_ij.shape[1]

    # メモリの関係上、分割して計算した方が早い。amsr_ij->amsr_ij2, modis_ij->modis_ij2 (分割前->分割後)
    length = imga_ij.shape[1]
    range_list = np.linspace(0, length, n_proc).astype(np.int32)

    for j in range(range_list.shape[0]-1):
        findex = range_list[j]
        lindex = range_list[j+1]
        imga_ij2 = imga_ij[:, findex:lindex]
        imgb_ij02 = imgb_ij[:, findex:lindex]

        imgb_ij2 = imgb_ij02+index_ij.T
        imgb_ij2 = imgb_ij2.transpose(2, 0, 1)
        n_imgb, n_kernel, _ = imgb_ij2.shape
        imgb_ij2 = imgb_ij2.reshape(n_imgb*n_kernel, 2)

        # modis_index2が配列ドメイン内という条件で、modis_index2上のマスクインデックスを取得
        mask_index = np.where((imgb_ij2[:, 0] > 0) & (imgb_ij2[:, 0] < imgb_nx) & (
            imgb_ij2[:, 1] > 0) & (imgb_ij2[:, 1] < imgb_ny))

        # マスクインデックス、modis_index,amsr_indexを使って各緯度経度値を取得
        imgb_ij3 = imgb_ij2[mask_index]
        lat_imgb = latb[imgb_ij3[:, 0], imgb_ij3[:, 1]]
        lon_imgb = lonb[imgb_ij3[:, 0], imgb_ij3[:, 1]]

        lon_imga = lona[imga_ij2[0], imga_ij2[1]]
        lat_imga = lata[imga_ij2[0], imga_ij2[1]]
        sp_imga = sp[:, imga_ij2[0], imga_ij2[1]]
        lon_imga = np.tile(lon_imga, (range_x*range_y*4, 1)
                           ).T.reshape(-1)[mask_index]
        lat_imga = np.tile(lat_imga, (range_x*range_y*4, 1)
                           ).T.reshape(-1)[mask_index]
        x_imga = np.tile(sp_imga[0], (range_x*range_y*4, 1)
                         ).T.reshape(-1)[mask_index]
        y_imga = np.tile(sp_imga[1], (range_x*range_y*4, 1)
                         ).T.reshape(-1)[mask_index]
        z_imga = np.tile(sp_imga[2], (range_x*range_y*4, 1)
                         ).T.reshape(-1)[mask_index]
        sp_imga = np.array([x_imga, y_imga, z_imga]).T

        p0 = np.array(geo_info.transform_lla_to_ecef(lat_imga, lon_imga, 0)).T
        b = p0-sp_imga
        i, j, k = calc_boresight_basis_vectors(p0, sp_imga)

        p = np.array(geo_info.transform_lla_to_ecef(lat_imgb, lon_imgb, 0)).T
        # print(p.shape,i.shape,j.shape,k.shape,b.shape,p0.shape)
        az, el = calc_local_az_el_angle(i, j, k, b, p0, p)
        az = np.int32(az*180/np.pi*100+offset)
        el = np.int32(el*180/np.pi*100+offset)
        coef_array0 = ap[el, az]

        # 結果をn_modis*n_kernelの配列を作成し、算出された係数を代入。その後、3次元に戻す。
        coef_array = np.full(n_imgb*n_kernel, np.NaN)
        coef_array[mask_index] = coef_array0
        coef_array = coef_array.reshape(n_imgb, n_kernel)
    #        print(coef_array.shape,modis_ij2.shape)

        # maskされたmodis_indexを使用してistを抽出。同じように3次元に戻す。そして、ist_amsrをnan値を考慮して算出。その後、amsrの座標系に戻す。
        imgb2 = np.full(n_imgb*n_kernel, np.NaN)
        imgb2[mask_index] = imgb[imgb_ij3[:, 0], imgb_ij3[:, 1]]
        imgb2 = imgb2.reshape(n_imgb, n_kernel)
        imgb_a0 = np.nansum(imgb2*coef_array, axis=1) / \
            np.nansum(coef_array, axis=1)
        imgb_a[imga_ij2[0], imga_ij2[1]] = imgb_a0

    return imgb_a


def interpolation_89G_AB(latb, lonb, tbb, lata, lona, spa, epsg, range_x=3, range_y=3, res=12000):
    """
    For polar region input epsg is 3413:arctic, 3976:antarctic
    coord index calculated from the identified epsg is only used for finding a neighborhood pixel.
    """
    # -------------get antenna pattern--------
    # 89GHz
    int_ap = 0.01
    integ_time = 1.3
    ifov_angle = 0.15
    antenna_func = get_antenna_pattern_gaussian_beam
    radius = 5
    # antenna_func=kpmw.get_antenna_pattern_bessel_beam
    az = np.arange(-radius, radius+int_ap, int_ap)
    el = np.arange(-radius, radius+int_ap, int_ap)
    X, Y = np.meshgrid(az, el)
    func_args = [ifov_angle, ifov_angle, X, Y]
    ap = antenna_pattern_integration(integ_time, 0.1, antenna_func, func_args)
    offset = radius/int_ap

    # ---------projection-----------
    # set transform
    source = osr.SpatialReference()
    source.ImportFromEPSG(4326)
    target = osr.SpatialReference()
    target.ImportFromEPSG(epsg)  # 3411:arctic, 3412:antarctic
    proj = geo_map.Search(source, target, latb, lonb, res)
    imgb_ij, imga_ij = proj.search_index(source, lata, lona)
    tbb2 = interpolation(latb, lonb, tbb, imgb_ij, lata,
                         lona, spa, imga_ij, ap, offset, range_x, range_y)
    return tbb2


def destriping_amsr2_89GHz(tba, lata, lona, tbb, latb, lonb, spa2, spb2, lata2, lona2, latb2, lonb2,
                           iter_n=250, param1=0.05, m_stop=6, m_stop2=5, xcut=10, epsg=3976):
    '''
    lata2, lona2, spa2 correspond to the matrixes of y-directional intermediate position for lata, lona2, spa2.
    intermediate positions are derived from _get_AMSR2_intermediate_positions.
    '''
    # ----------destriping for B-horn------------------
    ua0, param2a = img_destripe.proc_bouali_destriping(
        tba, iter_n, param1, m_stop)
    ub0, param2b = img_destripe.proc_bouali_destriping(
        tbb, iter_n, param1, m_stop)
    # --------interpolation into A-horn grid--------------------------
    tba2 = interpolation_89G_AB(lata, lona, ua0, latb2, lonb2, spb2, epsg=epsg)
    tbb2 = interpolation_89G_AB(latb, lonb, ub0, lata2, lona2, spa2, epsg=epsg)
    # ----------destriping for A-horn------------------
    length = tba.shape[0]+tbb2.shape[0]
    width = tba.shape[1]
    tba3 = np.zeros([length, width])
    tba3[::2, :] = tba
    tba3[1::2, :] = tbb2
    tba3 = tba3[:, xcut:-xcut]
    ua0, _ = img_destripe.proc_bouali_destriping(tba3, iter_n, param1, m_stop2)
    ua0 = ua0[0::2, :]
    ua = np.zeros(tba.shape)
    ua[:, xcut:-xcut] = ua0

    length = tbb.shape[0]+tba2.shape[0]
    width = tbb.shape[1]
    tbb3 = np.zeros([length, width])
    tbb3[::2, :] = tbb
    tbb3[1::2, :] = tba2
    ub0, _ = img_destripe.proc_bouali_destriping(tbb3, iter_n, param1, m_stop2)
    ub = ub0[0::2, :]

    return ua, ub


def destriping2_amsr2_89GHz(tba, lata, lona, tbb, latb, lonb, spa2, spb2, lata2, lona2, latb2, lonb2,
                            iter_n=250, param1=0.05, m_stop2=5, xcut=10, epsg=3976):
    '''
    lata2, lona2, spa2 correspond to the matrixes of y-directional intermediate position for lata, lona2, spa2.
    intermediate positions are derived from _get_AMSR2_intermediate_positions.
    '''
    # ----------destriping for B-horn------------------
    ua0 = tba
    ub0 = tbb
    # --------interpolation into A-horn grid--------------------------
    tba2 = interpolation_89G_AB(lata, lona, ua0, latb2, lonb2, spb2, epsg=epsg)
    tbb2 = interpolation_89G_AB(latb, lonb, ub0, lata2, lona2, spa2, epsg=epsg)
    # ----------destriping for A-horn------------------
    length = tba.shape[0]+tbb2.shape[0]
    width = tba.shape[1]
    tba3 = np.zeros([length, width])
    tba3[::2, :] = tba
    tba3[1::2, :] = tbb2
    tba3 = tba3[:, xcut:-xcut]
    ua0, _ = img_destripe.proc_bouali_destriping(tba3, iter_n, param1, m_stop2)
    ua0 = ua0[0::2, :]
    ua = np.zeros(tba.shape)
    ua[:, xcut:-xcut] = ua0

    length = tbb.shape[0]+tba2.shape[0]
    width = tbb.shape[1]
    tbb3 = np.zeros([length, width])
    tbb3[::2, :] = tbb
    tbb3[1::2, :] = tba2
    ub0, _ = img_destripe.proc_bouali_destriping(tbb3, iter_n, param1, m_stop2)
    ub = ub0[0::2, :]

    return ua, ub


def destriping3_amsr2_89GHz(tba, lata, lona, tbb, latb, lonb, spa2, spb2, lata2, lona2, latb2, lonb2,
                            iter_n=250, param1=0.05, m_stop=3, m_stop2=3, gain=0.25, xcut=10, epsg=3976):
    '''
    lata2, lona2, spa2 correspond to the matrixes of y-directional intermediate position for lata, lona2, spa2.
    intermediate positions are derived from _get_AMSR2_intermediate_positions.
    '''
    # ----------destriping for B-horn------------------
    ua0, param2a = img_destripe.proc_bouali_destriping(
        tba, iter_n, param1, m_stop, gain)
    ub0, param2b = img_destripe.proc_bouali_destriping(
        tbb, iter_n, param1, m_stop, gain)
    # --------interpolation into A-horn grid--------------------------
    tba2 = interpolation_89G_AB(lata, lona, ua0, latb2, lonb2, spb2, epsg=epsg)
    tbb2 = interpolation_89G_AB(latb, lonb, ub0, lata2, lona2, spa2, epsg=epsg)
    # ----------destriping for A-horn------------------
    length = tba.shape[0]+tbb2.shape[0]
    width = tba.shape[1]
    tba3 = np.zeros([length, width])
    tba3[::2, :] = ua0
    tba3[1::2, :] = tbb2
    tba3 = tba3[:, xcut:-xcut]
    ua0, _ = img_destripe.proc_bouali_destriping(
        tba3, iter_n, param1, m_stop2, gain)
    ua0 = ua0[0::2, :]
    ua = np.zeros(tba.shape)
    ua[:, xcut:-xcut] = ua0

    length = tbb.shape[0]+tba2.shape[0]
    width = tbb.shape[1]
    tbb3 = np.zeros([length, width])
    tbb3[::2, :] = ub0
    tbb3[1::2, :] = tba2
    ub0, _ = img_destripe.proc_bouali_destriping(
        tbb3, iter_n, param1, m_stop2, gain)
    ub = ub0[0::2, :]

    return ua, ub


def _get_AMSR2_intermediate_positions(lat, lon, orb, region, orb_i, orb_p, ecl):
    dorb_i = 7.5

    if orb == "D":
        dasc_lon = 0.8+180+12.25
    elif orb == "A":
        if region == "S":
            dasc_lon = 0.8+25
        elif region == "N":
            dasc_lon = 0.8

    orb_i = float(orb_i[0:6])+dorb_i
    orb_p = float(orb_p[0:4])/60/24  # [days]
    ecl = float(ecl)+dasc_lon
    source = osr.SpatialReference()
    source.ImportFromEPSG(4326)
    target = osr.SpatialReference()
    target_proj4 = geo_info.get_space_oblique_mercator_proj4(
        inc_angle=orb_i, ps_rev=orb_p, asc_lon=ecl, false_e=0, false_n=0)
    target.ImportFromProj4(target_proj4)

    coord = np.array([lat.reshape(-1), lon.reshape(-1)]).T
    coord_transform = osr.CoordinateTransformation(source, target)
    coord_itransform = osr.CoordinateTransformation(target, source)
    coord2 = np.array(coord_transform.TransformPoints(coord))[:, 0:2]
    coord3_y = coord2[:, 0].reshape(lat.shape)
    coord3_x = coord2[:, 1].reshape(lat.shape)
    coord3_y = (coord3_y[0:-1]+coord3_y[1:])/2
    coord3_x = (coord3_x[0:-1]+coord3_x[1:])/2
    coord = np.array([coord3_y.reshape(-1), coord3_x.reshape(-1)]).T
    coord = np.array(coord_itransform.TransformPoints(coord))[:, 0:2]
    lat2 = coord[:, 0].reshape(coord3_x.shape)
    lon2 = coord[:, 1].reshape(coord3_x.shape)
    return lat2, lon2
