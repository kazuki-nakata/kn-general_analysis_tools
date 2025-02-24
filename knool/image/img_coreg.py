import numpy as np
from ..fortlib import coreg_tool


def calc_offset_by_similarity(
    img1, img2, iy, ix, sim_type=int(1), ndw=int(13), dst=int(1), nsw=int(11), sst=int(1), disdep=float(0.5)
):
    """
    img1 and img2: master and slave images
    ix: 1D array of x coordinate where offset is calculated
    iy: 1D array of y coordinate where offset is calculated
    sim_type: 1=>zncc 2=>ssd 3=>sad 4=>ncc 5=>bbs
    ndw and dst : window size and stride
    nsw and sst : search window size and stride
    disdep: parameter for calculating bbs
    """
    sim, dy, dx = coreg_tool.calc_offset_type1(
        img1, img2, iy, ix, sim_type, ndw, dst, nsw, sst, disdep)
    return sim, dy, dx


def calc_offset_by_lucas_kanade(img1, img2, iy, ix, ndw=int(13), dst=int(1)):
    """
    multi band image can be processed.
    img1 and img2: master and slave images
    ix: 1D array of x coordinate where offset is calculated
    iy: 1D array of y coordinate where offset is calculated
    sim_type: 1=>zncc 2=>ssd 3=>sad 4=>ncc 5=>bbs
    ndw and dst : window size and stride
    """
    dy, dx = coreg_tool.lucas_kanade(img1, img2, iy, ix, ndw, dst)
    return dy, dx

# Horn-Schunck法の実装


def horn_schunck(Ix, Iy, It, alpha=1.0, num_iterations=100):
    u = np.zeros(Ix.shape)
    v = np.zeros(Ix.shape)

    avg = np.array([[1/12, 1/6, 1/12], [1/6, 0, 1/6], [1/12, 1/6, 1/12]])
    # avg = np.array([[1/24, 1/12, 1/24], [1/12, 0.5, 1/12], [1/24, 1/12, 1/24]])

    for _ in range(num_iterations):
        # u_avg = (np.roll(u, 1, axis=0) + np.roll(u, -1, axis=0) +
        #          np.roll(u, 1, axis=1) + np.roll(u, -1, axis=1)) / 4.0
        # v_avg = (np.roll(v, 1, axis=0) + np.roll(v, -1, axis=0) +
        #          np.roll(v, 1, axis=1) + np.roll(v, -1, axis=1)) / 4.0
        # u_avg = signal.convolve2d(u, avg, boundary='symm', mode='same')
        # v_avg = signal.convolve2d(v, avg, boundary='symm', mode='same')
        u_avg = filter2(u, avg)
        v_avg = filter2(v, avg)

        P = (Ix * u_avg + Iy * v_avg + It) / (alpha**2 + Ix**2 + Iy**2)
        u = u_avg - Ix * P
        v = v_avg - Iy * P

    return u, v


# Horn-Schunck法の実装


def horn_schunck_two_variables(Ix1, Iy1, It1, Ix2, Iy2, It2, alpha=1.0, num_iterations=100):
    u = np.zeros(Ix1.shape)
    v = np.zeros(Iy1.shape)
    avg = np.array([[1/12, 1/6, 1/12], [1/6, 0, 1/6], [1/12, 1/6, 1/12]])

    # ガウス・ザイデル法を用いて反復計算
    for _ in range(num_iterations):
        # u_avg = (np.roll(u, 1, axis=0) + np.roll(u, -1, axis=0) +
        #          np.roll(u, 1, axis=1) + np.roll(u, -1, axis=1)) / 4.0
        # v_avg = (np.roll(v, 1, axis=0) + np.roll(v, -1, axis=0) +
        #          np.roll(v, 1, axis=1) + np.roll(v, -1, axis=1)) / 4.0
        u_avg = filter2(u, avg)
        v_avg = filter2(v, avg)
        a1 = Ix1**2 + Ix2**2 + 2*(alpha**2)
        a2 = b1 = Ix1*Iy1+Ix2*Iy2
        a3 = 2*(alpha**2)*u_avg-Ix1*It1-Ix2*It2
        b2 = Iy1**2 + Iy2**2 + 2*(alpha**2)
        b3 = 2*(alpha**2)*v_avg-Iy1*It1-Iy2*It2
        u = (a3*b2-a2*b3)/(a1*b2-a2*b1)
        v = (a1*b3-a3*b1)/(a1*b2-a2*b1)
    return u, v
