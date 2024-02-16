import numpy as np
import cv2
from scipy import ndimage
import math


def get_psnr(image1, image2):
    PSNR, diff_square = cv2.quality.QualityPSNR_compute(image1, image2)
    return PSNR, diff_square


def get_ssim(image1, image2):
    SSIM, diff_square = cv2.quality.QualitySSIM_compute(image1, image2)
    return SSIM, diff_square


def get_mse(image1, image2):
    mse, diff_square = cv2.quality.QualityMSE_compute(image1, image2)
    return mse, diff_square


def get_snr(image1, image2):
    mse, diff_square = cv2.quality.QualityMSE_compute(image1, image2)
    snr = 10 * np.log10(np.var(image1) / mse[0])
    return snr, diff_square


def get_vifp_mscale(ref, dist):
    # H. R. Sheikh and A. C. Bovik, "Image Information and Visual Quality", IEEE Transactions on Image Processing
    sigma_nsq = 1000
    eps = 1e-10

    num = 0.0
    den = 0.0
    for scale in range(1, 5):

        N = 2 ** (4 - scale + 1) + 1
        sd = N / 5.0

        if scale > 1:
            ref = ndimage.gaussian_filter(ref, sd)
            dist = ndimage.gaussian_filter(dist, sd)
            ref = ref[::2, ::2]
            dist = dist[::2, ::2]

        mu1 = ndimage.gaussian_filter(ref, sd)
        mu2 = ndimage.gaussian_filter(dist, sd)
        mu1_sq = mu1 * mu1
        mu2_sq = mu2 * mu2
        mu1_mu2 = mu1 * mu2
        sigma1_sq = ndimage.gaussian_filter(ref * ref, sd) - mu1_sq
        sigma2_sq = ndimage.gaussian_filter(dist * dist, sd) - mu2_sq
        sigma12 = ndimage.gaussian_filter(ref * dist, sd) - mu1_mu2

        sigma1_sq[sigma1_sq < 0] = 0
        sigma2_sq[sigma2_sq < 0] = 0

        g = sigma12 / (sigma1_sq + eps)
        sv_sq = sigma2_sq - g * sigma12

        g[sigma1_sq < eps] = 0
        sv_sq[sigma1_sq < eps] = sigma2_sq[sigma1_sq < eps]
        sigma1_sq[sigma1_sq < eps] = 0

        g[sigma2_sq < eps] = 0
        sv_sq[sigma2_sq < eps] = 0

        sv_sq[g < 0] = sigma2_sq[g < 0]
        g[g < 0] = 0
        sv_sq[sv_sq <= eps] = eps

        num += np.sum(np.log10(1 + g * g * sigma1_sq / (sv_sq + sigma_nsq)))
        den += np.sum(np.log10(1 + sigma1_sq / sigma_nsq))

    try:
        vifp = math.log10(num / den)
    except:
        vifp = 999

    return vifp


def get_directional_total_spectrum_ratio(ori_img, dist_img, sampling_freq, axis=0):
    """
    in the case of all the freq sampling, you can input sampling_freq = np.fft.fftfreq(ori_img.shape[axis], d=1)[:int(ori_img.shape[axis]/2)]
    """
    num = ori_img.shape[axis]
    if axis == 1:
        ori_img = ori_img.T
        dist_img = dist_img.T

    ori_img2 = ori_img-np.mean(ori_img, axis=0)
    dist_img2 = dist_img-np.mean(dist_img, axis=0)

    # 0:y directional spectrums
    ori_pow = np.abs(np.fft.fft(ori_img2, axis=0)[:int(num/2+1)])**2
    dist_pow = np.abs(np.fft.fft(dist_img2, axis=0)[
                      :int(num/2+1)])**2  # 0:y directional spectrums

    ori_pow = np.mean(ori_pow, axis=1)
    dist_pow = np.mean(dist_pow, axis=1)
    sampling_freq_num = (sampling_freq*num).astype(np.int32)
    ori_pow2 = ori_pow[sampling_freq_num]
    dist_pow2 = dist_pow[sampling_freq_num]
    ratio = np.sum(dist_pow2)/np.sum(ori_pow2)
    return ratio


def get_directional_distorsion_index(ori_img, dist_img, sampling_freq, axis=0):
    """
    in the case of all the freq sampling, you can input sampling_freq = np.fft.fftfreq(ori_img.shape[axis], d=1)[:int(ori_img.shape[axis]/2)]
    """
    num = ori_img.shape[axis]
    if axis == 1:
        ori_img = ori_img.T
        dist_img = dist_img.T

    ori_img2 = ori_img-np.mean(ori_img, axis=0)
    dist_img2 = dist_img-np.mean(dist_img, axis=0)

    # 0:y directional spectrums
    ori_pow = np.abs(np.fft.fft(ori_img2, axis=0)[:int(num/2+1)])**2
    dist_pow = np.abs(np.fft.fft(dist_img2, axis=0)[
                      :int(num/2+1)])**2  # 0:y directional spectrums

    ori_pow = np.mean(ori_pow, axis=1)
    dist_pow = np.mean(dist_pow, axis=1)
    sampling_freq_num = (sampling_freq*num).astype(np.int32)
    ori_pow2 = ori_pow[sampling_freq_num]
    dist_pow2 = dist_pow[sampling_freq_num]
    di = 1-np.sum(np.abs(ori_pow2-dist_pow2)/ori_pow2)/sampling_freq.shape[0]
    return di
