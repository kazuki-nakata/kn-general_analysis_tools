import numpy as np

def image_split(self,img, n, m, b):
    print('split a image')
    print(img.shape[1])
    print(img.shape[2])
    pv_size = img.shape[1] // m
    ph_size = img.shape[2] // n
    v_size = img.shape[1] // pv_size * pv_size
    h_size = img.shape[2] // ph_size * ph_size
    img = img[0:b, :v_size, :h_size]
    out_img = []
    [out_img.extend(np.split(h_img, n, axis=2))
        for h_img in np.split(img, m, axis=1)]
    return out_img

def image_reunion(self,img, n, m, ijksize):
    out_img = np.zeros(ijksize, dtype=img[0].dtype)
    h = img[0].shape[2]
    v = img[0].shape[1]
    print(h)
    print(v)
    for vnum in range(m):
        for hnum in range(n):
            n1 = vnum * n + hnum
            hm = h * hnum
            hl = h * (hnum + 1)
            vm = v * vnum
            vl = v * (vnum + 1)
            out_img[:, vm:vl, hm:hl] = img[n1][:, :, :]
    return out_img