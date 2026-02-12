import numpy as np
from ..fortlib import deblur
from scipy import signal
from ..satellite.algorithm import pmw_processor
from scipy.signal import convolve2d
from scipy.sparse import lil_matrix


def Tikhonov_L2(G, alpha):
    U, S, V = np.linalg.svd(G, full_matrices=True)
    Sigma = np.zeros((G.shape[0], G.shape[1]))
    square_len = min((G.shape[0], G.shape[1]))
    Sigma[:square_len, :square_len] = np.diag(S/(S*S+alpha))
    Gg = V.T.dot(Sigma.dot(U.T))
    return Gg


def TruncatedSVD(G, n_elements):
    U, S, V = np.linalg.svd(G, full_matrices=True)
    Sigma = np.zeros((G.shape[0], G.shape[1]))
    square_len = min((G.shape[0], G.shape[1]))
    Sigma[:square_len, :square_len] = np.diag(1/S)
    U = U[:, :n_elements]
    Sigma = Sigma[:n_elements, :n_elements]
    V = V[:n_elements, :]
    Gg = V.T.dot(Sigma.dot(U.T))
    return Gg


def backus_gilbert(ap_s, ap_t, noise, gamma, interval=1, mode="SOLA"):
    """
    SOLA: Subtractive Optimally Localized Averages
    MOLA: Multiplicative Optimally Localized Averages
    """

    if mode == "SOLA":
        size_covm = np.eye(ap_s.shape[0], ap_s.shape[0])*noise
        G = ap_s.dot(ap_s.T)*interval*interval
        u = np.sum(ap_s, axis=1)*interval*interval
        v = ap_s.dot(ap_t)*interval*interval
        Z = G*np.cos(gamma)+size_covm*np.sin(gamma)
        Zinv = np.linalg.inv(Z)
        beta = (1-u.T.dot(Zinv.dot(v)*np.cos(gamma)))/(u.T.dot(Zinv.dot(u)))
        filter = Zinv.dot(v*np.cos(gamma)+beta*u)

    return filter


def run_rsir(grid_x, grid_y, mask, ap3, wsize, res, fwhm, w, val, iterate):
    output = deblur.sir2(grid_x, grid_y, val, mask, ap3,
                         wsize, res, fwhm, w, iterate)
    return output.T


def run_mart(grid_x, grid_y, mask, ap3, wsize, res, fwhm, w, frmin, val, iterate):
    output = deblur.mart(grid_x, grid_y, val, mask, ap3,
                         wsize, res, fwhm, w, iterate, frmin)
    return output.T


def run_bgrad(grid_x, grid_y, ap3, wsize, res, fwhm, w, w2, val, iterate):
    output = deblur.banach_gradient(
        grid_x, grid_y, val, ap3, wsize, res, fwhm, w, w2, iterate)
    return output.T


def convert_4DFilter_to_4DInverseFilter(aarray0, nh):
    ap2 = np.zeros([nh*2+1, nh*2+1, nh*2+1, nh*2+1])
    for i in range(-nh, nh+1, 1):
        for j in range(-nh, nh+1, 1):
            if j > 0:
                j_t = (j, 0)
            else:
                j_t = (0, -j)
            if i > 0:
                i_t = (i, 0)
            else:
                i_t = (0, -i)
            tmp = np.pad(aarray0[j+nh, i+nh, :, :],
                         pad_width=(j_t, i_t), mode="edge")
            length, width = tmp.shape
            if j > 0:
                jmin = 0
                jmax = length-j
            else:
                jmin = -j
                jmax = length
            if i > 0:
                imin = 0
                imax = width-i
            else:
                imin = -i
                imax = width

            tmp2 = tmp[jmin:jmax, imin:imax]
            height2, width2 = tmp.shape
            tmp2[0:jmin, :] = tmp2[0:jmin, :]+tmp[0:jmin, imin:imax][::-1, :]
            tmp2[:, 0:imin] = tmp2[:, 0:imin]+tmp[jmin:jmax, 0:imin][:, ::-1]
            if jmax != height2:
                tmp2[-(height2-jmax):, :] = tmp2[-(height2-jmax):, :] + \
                    tmp[jmax:height2, imin:imax][::-1, :]
            if imax != width2:
                tmp2[:, -(width2-imax):] = tmp2[:, -(width2-imax):] + \
                    tmp[jmin:jmax, imax:width2][:, ::-1]
            tmp2 = tmp2/np.sum(tmp2)
            ap2[j+nh, i+nh, :, :] = tmp2

    return ap2


def create_convolution_matrix_from_2DFilter(array0, ratio=1, padding=0, reg=True):
    '''output shape is (point_y,point_x,filter_y,filter_x)'''
    nh = int((array0.shape[0]-1)/2)
    nh2 = (nh-padding)*ratio
    if nh2 % int(nh2) == 0:
        nh2 = int(nh2)
        ap2 = np.zeros([nh2*2+1, nh2*2+1, nh*2+1, nh*2+1])
        print(ap2.shape,-nh+padding,nh+1-padding)
        for i in range(-nh+padding, nh+1-padding, int(1/ratio)):
            for j in range(-nh+padding, nh+1-padding, int(1/ratio)):
                if j > 0:
                    j_t = (j, 0)
                else:
                    j_t = (0, -j)
                if i > 0:
                    i_t = (i, 0)
                else:
                    i_t = (0, -i)
                tmp = np.pad(array0[:, :],
                             pad_width=(j_t, i_t), mode="edge")
                length, width = tmp.shape
                if j > 0:
                    jmin = 0
                    jmax = length-j
                else:
                    jmin = -j
                    jmax = length
                if i > 0:
                    imin = 0
                    imax = width-i
                else:
                    imin = -i
                    imax = width

                tmp2 = tmp[jmin:jmax, imin:imax]
                height2, width2 = tmp.shape
                tmp2[0:jmin, :] = tmp2[0:jmin, :] + \
                    tmp[0:jmin, imin:imax][::-1, :]
                tmp2[:, 0:imin] = tmp2[:, 0:imin] + \
                    tmp[jmin:jmax, 0:imin][:, ::-1]
                if jmax != height2:
                    tmp2[-(height2-jmax):, :] = tmp2[-(height2-jmax):, :] + \
                        tmp[jmax:height2, imin:imax][::-1, :]
                if imax != width2:
                    tmp2[:, -(width2-imax):] = tmp2[:, -(width2-imax):] + \
                        tmp[jmin:jmax, imax:width2][:, ::-1]
                if reg:
                    tmp2 = tmp2/np.sum(tmp2)
                # print(tmp2.shape,int(j/ratio)+nh2,int(i/ratio)+nh2)
                ap2[int(j*ratio)+nh2, int(i*ratio)+nh2, :, :] = tmp2
    else:
        print("ratio should be ajusted.")
    return ap2.reshape((nh2*2+1)**2, (nh*2+1)**2)


def convert_2DFilter_to_4DInverseFilter(aarray0, nh):
    '''remove!!'''
    ap2 = np.zeros([nh*2+1, nh*2+1, nh*2+1, nh*2+1])
    for i in range(-nh, nh+1, 1):
        for j in range(-nh, nh+1, 1):
            if j > 0:
                j_t = (j, 0)
            else:
                j_t = (0, -j)
            if i > 0:
                i_t = (i, 0)
            else:
                i_t = (0, -i)
            tmp = np.pad(aarray0[:, :],
                         pad_width=(j_t, i_t), mode="edge")
            length, width = tmp.shape
            if j > 0:
                jmin = 0
                jmax = length-j
            else:
                jmin = -j
                jmax = length
            if i > 0:
                imin = 0
                imax = width-i
            else:
                imin = -i
                imax = width

            tmp2 = tmp[jmin:jmax, imin:imax]
            height2, width2 = tmp.shape
            tmp2[0:jmin, :] = tmp2[0:jmin, :]+tmp[0:jmin, imin:imax][::-1, :]
            tmp2[:, 0:imin] = tmp2[:, 0:imin]+tmp[jmin:jmax, 0:imin][:, ::-1]
            if jmax != height2:
                tmp2[-(height2-jmax):, :] = tmp2[-(height2-jmax):, :] + \
                    tmp[jmax:height2, imin:imax][::-1, :]
            if imax != width2:
                tmp2[:, -(width2-imax):] = tmp2[:, -(width2-imax):] + \
                    tmp[jmin:jmax, imax:width2][:, ::-1]
            tmp2 = tmp2/np.sum(tmp2)
            ap2[j+nh, i+nh, :, :] = tmp2

    return ap2


def convert_3DFilter_to_4DInverseFilter(test, nh, point):
    nh = 20
    nh2 = nh*2
    line = 20
    nsize = nh*2+1
    for inum in range(point-nh, point-nh+1):
        print(inum)
        ap2 = np.zeros([nh+1, nh+1, nh+1, nh+1])
        for j in range(41):
            line2 = line-nh+j
            for i in range(nh2):
                inum2 = inum-nh+i
                test2 = np.zeros([nsize, nsize])
                jmin = j-nh
                jmax = j+nh
                imin = i-nh
                imax = i+nh
                if jmax > nh2:
                    jmax = nh2
                if jmin < 0:
                    jmin = 0
                if imax > nh2:
                    imax = nh2
                if imin < 0:
                    imin = 0
                jmin2 = nh-j
                jmax2 = nh*2-j
                imin2 = nh-i
                imax2 = nh+nh-i
                if jmax2 > nh2:
                    jmax2 = nh2
                if jmin2 < 0:
                    jmin2 = 0
                if imax == nh2:
                    imin2 = 0
                    imax2 = nh*2-imin
                if imin == 0:
                    imin2 = nh*2-imax
                    imax2 = nh2
                if (imin == 0) & (imax == nh2):
                    imin2 = 0
                    imax2 = nh2
                if jmax == nh2:
                    jmin2 = 0
                    jmax2 = nh*2-jmin
                if jmin == 0:
                    jmin2 = nh*2-jmax
                    jmax2 = nh2
                if (jmin == 0) & (jmax == 40):
                    jmin2 = 0
                    jmax2 = nh2

                test2[jmin2:jmax2, imin2:imax2] = test[jmin:jmax, imin:imax, inum2]
                ap2[:, :, j, i] = test2
    return ap2


def get_ideal_antenna_pattern(nh_x, nh_y, fwhm_x, fwhm_y):
    antenna_func = pmw_processor.get_antenna_pattern_gaussian_beam
    az = np.arange(-nh_x, nh_x+1, 1)
    el = np.arange(-nh_y, nh_y+1, 1)
    X, Y = np.meshgrid(az, el)
    func_args = [fwhm_x, fwhm_y, X, Y]
    ap = antenna_func(*func_args)
    ap = ap/np.sum(ap)
    return ap


def estimate_resolution(D, smin, smax, sint, metric="rmse"):
    """D=GgG matrix,metric=rmse or corr"""
    height, width = D.shape
    corr_max = 0
    rmse_opt = 9.9E33
    for fwhm_y in np.arange(smin, smax+sint, sint):
        for fwhm_x in np.arange(smin, smax+sint, sint):
            G = get_ideal_antenna_pattern(
                (width-1)/2, (height-1)/2, fwhm_x, fwhm_y)
            # conv = signal.convolve2d(F, G, boundary='symm', mode='same')
            if metric == "rmse":
                diff = G.reshape(-1)-D.reshape(-1)
                rmse = np.sqrt(np.sum(diff**2)/height/width)
                if rmse < rmse_opt:
                    fwhm_opt_x = fwhm_x
                    fwhm_opt_y = fwhm_y
                    rmse_opt = rmse
                    gt = G
            if metric == "corr":
                corr = np.corrcoef(G.reshape(-1), D.reshape(-1))[0, 1]
                if corr > corr_max:
                    fwhm_opt_x = fwhm_x
                    fwhm_opt_y = fwhm_y
                    corr_max = corr
                    gt = G

    if metric == "rmse":
        metric_val = rmse_opt
    else:
        metric_val = corr_max

    return fwhm_opt_x, fwhm_opt_y, metric_val, gt


def create_convolution_matrix_sparse(G, nh, ratio=1):
    N = int(nh*2*ratio+1)
    N2 = nh*2+1
    A = lil_matrix((N*N, N2*N2))  # 大きさは [H*W, H*W]
    # 各 unit vector に対する畳み込み結果を A の各行として格納
    j = 0
    for i in range(0, N2*N2, int(1/ratio)):
        unit_img = np.zeros((N2, N2))
        unit_img[np.unravel_index(i, (N2, N2))] = 1.0
        conv_result = convolve2d(unit_img, G, mode='same', boundary='wrap')
        A[j, :] = conv_result.flatten()
        j = j+1
    return A.tocsr()
