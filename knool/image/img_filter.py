import numpy as np
import cv2
from scipy import ndimage


def gaussian_filter(img, sigma, mode="cv2"):
    if mode == "cv2":
        kernel_size = 6 * sigma - 1
        img2 = cv2.GaussianBlur(img, (kernel_size, kernel_size), sigma)
    elif mode == "scipy":
        ndimage.gaussian_filter(img, sigma)
    return img2
