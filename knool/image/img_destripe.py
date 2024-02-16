import numpy as np


def bouali_destriping(img, iter_n, param1=0.5, param2="TVbase", eps1=0.01, eps2=0.01):
    """# Bouali (2011) MODIS destriping, pattern 1"""
    grad_oa_x = np.zeros(img.shape)
    grad_oa_y = np.zeros(img.shape)
    grad_ua_x = np.zeros(img.shape)
    grad_ua_y = np.zeros(img.shape)
    add1 = np.zeros(img.shape)
    add2 = np.zeros(img.shape)

    ua = img
    grad_oa_x[:, :-1] = img[:, 1:]-img[:, :-1]
    grad_oa_y[:-1, :] = img[1:, :]-img[:-1, :]

    if param2 == "TVbase":
        param2 = np.mean(np.abs(grad_oa_y))/np.mean(np.abs(grad_oa_x))

    for i in range(iter_n):
        grad_ua_x[:, :-1] = ua[:, 1:]-ua[:, :-1]
        g1 = grad_ua_x-grad_oa_x
        g1 = g1/np.sqrt(g1**2+eps1)
        add1[:, 1:] = g1[:, 1:]-g1[:, :-1]

        grad_ua_y[:-1, :] = ua[1:, :]-ua[:-1, :]
        g1 = grad_ua_y/np.sqrt(grad_ua_y**2+eps2)
        add2[1:, :] = g1[1:, :]-g1[:-1, :]
        ua = ua+param1*(add1+param2*add2)

    return ua


def proc_bouali_destriping(img, iter_n, param1=0.5, n_stop=10, param2="Auto", gain=1,
                           img_quality_method=None, iqm_param=None, thres=0.95, eps1=0.01, eps2=0.01):
    """
    Cascade iteration by Bouali (2011) MODIS destriping.
    iter_n and param1 is iteration parameters.
    base of image quality and its threshold is get_directional_distorsion_index and 0.95, respectively.
    initial parameter for y-directional TV (param2) is automatically determined by param2=np.mean(np.abs(grad_oa_y))/np.mean(np.abs(grad_oa_x))
    In the case of use of img_quality_method (iqm), you must input iqm_param (arg of iqm) and threshold.
    """
    uk = np.zeros(img.shape)
    grad_oa_x = np.zeros(img.shape)
    grad_oa_y = np.zeros(img.shape)
    grad_oa_x[:, :-1] = img[:, 1:]-img[:, :-1]
    grad_oa_y[:-1, :] = img[1:, :]-img[:-1, :]
    i = 0
    qi = 0
    if param2 == "Auto":
        param2 = np.mean(np.abs(grad_oa_y))/np.mean(np.abs(grad_oa_x)) * gain

    param2_init = param2

    while qi < thres:
        i = i+1
        img2 = img-uk
        param2a = param2/(2**i)
        ua = bouali_destriping(img2, iter_n, param1=param1,
                               param2=param2a, eps1=eps1, eps2=eps2)
        uk = uk+ua
        if img_quality_method is not None:
            qi = img_quality_method(img, uk, *iqm_param)
        if i == n_stop:
            break
    return uk, param2_init


def tikhonov_L2_destripe(img, param1=0.5):  # direct Tikhnov L2 destriping
    img2 = np.mean(img, axis=1)
    dcomp = np.zeros(img2.shape)
    dcomp[0] = -1
    dcomp[1] = 1
    D = np.array([np.roll(dcomp, i)for i in range(dcomp.shape[0]-1)])
    # D[-1,:]=0
    print(D.shape)
    g1 = np.linalg.inv(np.dot(D.T, D)+param1*np.eye(img2.shape[0]))
    print(D.shape)
    g2 = np.dot(D.T, D).dot(img2)
    b = np.dot(g1, g2)
    print(b.shape, g1.shape, g2.shape)
    uk = img-b.reshape([b.shape[0], 1])
    return uk


def _svd_Tikhonov_L2_for_destripe(ap, alpha):
    U, S, V = np.linalg.svd(ap, full_matrices=True)
    Sigma = np.zeros((ap.shape[1], ap.shape[1]))
    square_len = min((ap.shape[0], ap.shape[1]))
    Sigma[:square_len, :square_len] = np.diag(S*S/(S*S+alpha))
    Gg = V.T.dot(Sigma).dot(V)
    return Gg


# Tikhnov L2 destriping based on SVD
def tikhonov_L2_destriping_SVD(img, param1=0.5):
    img2 = np.mean(img, axis=1)
    dcomp = np.zeros(img2.shape)
    dcomp[0] = -1
    dcomp[1] = 1
    D = np.array([np.roll(dcomp, i)for i in range(dcomp.shape[0]-1)])
    Gg = _svd_Tikhonov_L2_for_destripe(D, param1)
    print(Gg.shape)
    b = np.dot(Gg.T, img2)
    uk = img-b.reshape([b.shape[0], 1])
    return uk
