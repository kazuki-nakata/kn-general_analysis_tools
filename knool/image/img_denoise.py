import numpy as np
import pywt
import matplotlib.pyplot as plt
from skimage.restoration import denoise_tv_chambolle, denoise_tv_bregman


def image_normalization(src_img):
    """
    白飛び防止のための正規化処理
    cv2.imshowでwavelet変換された画像を表示するときに必要（大きい値を持つ画像の時だけ）
    """
    norm_img = (src_img - np.min(src_img)) / \
        (np.max(src_img) - np.min(src_img))
    return norm_img


def merge_images(cA, cH_V_D):
    """numpy.array を４つ(左上、(右上、左下、右下))連結させる"""
    cH, cV, cD = cH_V_D
    cH = image_normalization(cH)  # 外してもok
    cV = image_normalization(cV)  # 外してもok
    cD = image_normalization(cD)  # 外してもok
    # 元画像が2の累乗でない場合、端数ができることがあるので、サイズを合わせる。小さい方に合わせます。
    cA = cA[0:cH.shape[0], 0:cV.shape[1]]
    # 左上、右上、左下、右下、で画素をくっつける
    return np.vstack((np.hstack((cA, cH)), np.hstack((cV, cD))))


def coeffs_visualization(cof):
    norm_cof0 = cof[0]
    norm_cof0 = image_normalization(norm_cof0)  # 外してもok
    merge = norm_cof0
    for i in range(1, len(cof)):
        merge = merge_images(merge, cof[i])  # ４つの画像を合わせていく
    plt.figure(figsize=(10, 10))
    plt.imshow(merge, cmap="seismic")
    plt.colorbar()
    plt.show()


def func_soft_threshold(coeff, threshold):
    return np.sign(coeff) * np.maximum(np.abs(coeff) - threshold, 0)


def func_hard_threshold(coeff, threshold):
    return coeff * (np.abs(coeff) > threshold)


def thresholding(coeffs, threshold, threshold_func, level=2):
    cA, *detail_coeffs = coeffs
    detail_coeffs_corr = [(threshold_func(cH, threshold), threshold_func(cV, threshold), threshold_func(cD, threshold))
                          for cH, cV, cD in detail_coeffs[-level:]]
    print(len(coeffs), len(detail_coeffs), len(
        detail_coeffs_corr), len(detail_coeffs[:-level]))
    return [cA] + detail_coeffs[:-level] + detail_coeffs_corr


def wv_sum(coeff1, coeff2, gain=1):
    output = []
    cA1, *detail_coeffs1 = coeff1
    cA2, *detail_coeffs2 = coeff2
    output = [cA1+gain*cA2]
    for i in range(len(detail_coeffs1)):
        output.append((detail_coeffs1[i][0]+gain*detail_coeffs2[i][0], detail_coeffs1[i]
                      [1]+gain*detail_coeffs2[i][1], detail_coeffs1[i][2]+gain*detail_coeffs2[i][2]))
    return output


def ISTA(y, lam=50e-1, n_iter=200):
    level = 2
    wavelet = 'db1'  # Daubechies wavelet
    coeffs = pywt.wavedec2(y, wavelet, level=level)
    for k in range(n_iter):
        # convolve2d(np.sum(x, axis=0), ap, mode='same', boundary='symm')
        y_est = pywt.waverec2(coeffs, wavelet)
        bp = y_est-y  # convolve2d(fp - b, ap, mode='same', boundary='symm')
        coeffs = wv_sum(coeffs, pywt.wavedec2(
            bp, wavelet, level=level), gain=-1)
        cA, *detail_coeffs = coeffs
        coeffs = [cA] + [(func_soft_threshold(cH, lam), func_soft_threshold(cV, lam),
                          func_soft_threshold(cD, lam)) for cH, cV, cD in detail_coeffs]
    return pywt.waverec2(coeffs, wavelet)


def tv_bregman(y, weight=0.5, max_num_iter=5000, isotropic=False):
    return denoise_tv_bregman(y, weight=weight, max_num_iter=max_num_iter, isotropic=isotropic)
