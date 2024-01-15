import numpy as np
import multiprocessing as mp
from multiprocessing.connection import wait
import numpy as np
import cv2
import math
from PIL import Image
from ..fortlib import coreg_tool


def calc_offset(
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


def resize(img, pixel_x, pixel_y, resample=Image.BOX, box=None, reducing_gap=None):
    #     Image.NEAREST Image.BOX Image.BILINEAR Image.HAMMING Image.BICUBIC Image.LANCZOS
    img2 = Image.fromarray(img)
    img2 = img2.resize((pixel_y, pixel_x), resample=resample,
                       box=None, reducing_gap=None)
    return np.array(img2)


def split(img, n, m, b):
    print("split a image")
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


def reunion(img, n, m, ijksize):
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


def perspective_transformation(img, rot_x, rot_y, rot_z, meter, move):

    length = img.shape[1]
    width = img.shape[0]
    R = 5.0 / (meter / width)

    # ul = [0,0,0,1]
    # ur = [0,width,0,1]
    # bl = [length,0,0,1]
    # br = [length,width,0,1]

    ul = [-length / 2, -width / 2, 0, 1]
    ur = [-length / 2, width / 2, 0, 1]
    bl = [length / 2, -width / 2, 0, 1]
    br = [length / 2, width / 2, 0, 1]

    base = [0, 0, 0, 1]
    move2 = [0, 0, move, 1]

    rot_x = math.radians(rot_x)
    rot_y = math.radians(rot_y)
    rot_z = math.radians(rot_z)

    # X = np.array([[1,0,0],[0,math.cos(rot_x),math.sin(rot_x)],[0,-math.sin(rot_x),math.cos(rot_x)]],dtype=np.float32)
    # Y = np.array([[math.cos(rot_y),0,-math.sin(rot_y)],[0,1,0],[math.sin(rot_y),0,math.cos(rot_y)]],dtype=np.float32)
    # Z = np.array([[math.cos(rot_z),math.sin(rot_z),0],[-math.sin(rot_z),math.cos(rot_z),0],[0,0,1]],dtype=np.float32)

    X = np.array(
        [
            [1, 0, 0, 0],
            [0, math.cos(rot_x), math.sin(rot_x), 0],
            [0, -math.sin(rot_x), math.cos(rot_x), 0],
            [0, 0, 0, 1],
        ],
        dtype=np.float32,
    )
    Y = np.array(
        [
            [math.cos(rot_y), 0, -math.sin(rot_y), 0],
            [0, 1, 0, 0],
            [math.sin(rot_y), 0, math.cos(rot_y), 0],
            [0, 0, 0, 1],
        ],
        dtype=np.float32,
    )
    Z = np.array(
        [
            [math.cos(rot_z), math.sin(rot_z), 0, 0],
            [-math.sin(rot_z), math.cos(rot_z), 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1],
        ],
        dtype=np.float32,
    )

    ul_t = np.dot(np.dot(X, np.dot(Y, Z)), ul)[
        0:2] / (1 - np.dot(np.dot(X, np.dot(Y, Z)), ul)[2] / R)
    ur_t = np.dot(np.dot(X, np.dot(Y, Z)), ur)[
        0:2] / (1 - np.dot(np.dot(X, np.dot(Y, Z)), ur)[2] / R)
    bl_t = np.dot(np.dot(X, np.dot(Y, Z)), bl)[
        0:2] / (1 - np.dot(np.dot(X, np.dot(Y, Z)), bl)[2] / R)
    br_t = np.dot(np.dot(X, np.dot(Y, Z)), br)[
        0:2] / (1 - np.dot(np.dot(X, np.dot(Y, Z)), br)[2] / R)
    base_t = np.dot(np.dot(X, np.dot(Y, Z)), base)[
        0:2] / (1 - np.dot(np.dot(X, np.dot(Y, Z)), base)[2] / R)
    move_t = np.dot(np.dot(X, np.dot(Y, Z)), move2)[
        0:2] / (1 - np.dot(np.dot(X, np.dot(Y, Z)), move2)[2] / R)

    xmin = min(ul_t[1], ur_t[1], bl_t[1], br_t[1])
    ymin = min(ul_t[0], ur_t[0], bl_t[0], br_t[0])
    xmax = max(ul_t[1], ur_t[1], bl_t[1], br_t[1]) - xmin
    ymax = max(ul_t[0], ur_t[0], bl_t[0], br_t[0]) - ymin

    offset = np.array([-ymin, -xmin], dtype=np.float32)
    ul_t = ul_t + offset
    ur_t = ur_t + offset
    bl_t = bl_t + offset
    br_t = br_t + offset
    base_t = base_t + offset
    move_t = move_t + offset

    #    print(base_t)
    #    print(move_t)
    area = cv2.contourArea(np.float32([ul_t, ur_t, br_t, bl_t]))

    pts1 = np.float32([[0, 0], [0, width], [length, 0], [length, width]])
    pts2 = np.float32([ul_t, ur_t, bl_t, br_t])
    M = cv2.getPerspectiveTransform(pts1, pts2)

    #    M = np.delete(np.delete(np.dot(X,np.dot(Y,Z)),2,0),2,1)
    #    M = np.dot(X,np.dot(Y,Z))
    #    print(M)

    cols = int(ymax)
    rows = int(xmax)

    trans_img = cv2.warpPerspective(
        img, M, (cols, rows), borderMode=cv2.BORDER_REPLICATE)  # [::-1,:]
    #    trans_img = cv2.warpPerspective(img,M,(cols,rows))
    #    trans_img = cv2.getRectSubPix(trans_img,(400,400),(0,0))
    #    print("test",img.shape,trans_img.shape)
    #    print(trans_img[int(base_t[0]) - 10:int(base_t[0]) + 10, int(base_t[1]) - 10:int(base_t[1]) + 10, 0])
    # ----------------confirm move vector---------------------------------------------------------------------
    #     trans_img[int(base_t[1]) - 10:int(base_t[1]) + 10, int(base_t[0]) - 10:int(base_t[0]) + 10, 0] = 255
    #     trans_img[int(base_t[1]) - 10:int(base_t[1]) + 10, int(base_t[0]) - 10:int(base_t[0]) + 10, 1] = 0
    #     trans_img[int(base_t[1]) - 10:int(base_t[1]) + 10, int(base_t[0]) - 10:int(base_t[0]) + 10, 2] = 0
    #     trans_img[int(move_t[1]) - 10:int(move_t[1]) + 10, int(move_t[0]) - 10:int(move_t[0]) + 10, 0] = 0
    #     trans_img[int(move_t[1]) - 10:int(move_t[1]) + 10, int(move_t[0]) - 10:int(move_t[0]) + 10, 1] = 0
    #     trans_img[int(move_t[1]) - 10:int(move_t[1]) + 10, int(move_t[0]) - 10:int(move_t[0]) + 10, 2] = 255
    # -----------------------------------------------------------------------------------------------------------
    #    print(trans_img[int(base_t[0])-10:int(base_t[0])+10,int(base_t[1])-10:int(base_t[1])+10,0])
    return trans_img, area, pts1, pts2, base_t, move_t


def perspective_transformation_inverse(img, pts1, pts2):
    rows = int(max(pts2[0][1], pts2[1][1], pts2[2][1], pts2[3][1]))
    cols = int(max(pts2[0][0], pts2[1][0], pts2[2][0], pts2[3][0]))

    M = cv2.getPerspectiveTransform(pts1, pts2)
    trans_img = cv2.warpPerspective(
        img, M, (cols, rows), borderMode=cv2.BORDER_REPLICATE)
    return trans_img


def motion_blur(img, base, motion):
    dx = motion[0] - base[0]
    dy = motion[1] - base[1]
    theta = -math.atan(dy / dx) * 180 / 3.14
    alpha = dy / dx

    width = img.shape[1]
    length = img.shape[0]
    band = img.shape[2]
    size = int(math.sqrt(dx**2 + dy**2))

    if theta > 0:
        #       blur_kernel=np.full(size,1/float(size),dtype=np.float32)
        blur_kernel_y = np.zeros(size, dtype=np.int)
        blur_kernel_x = np.arange(size)
        num = 0
        for x in range(size):
            blur_kernel_y[num] = int(alpha * x)
            num = 1 + num
        ysize = np.max(blur_kernel_y) - np.min(blur_kernel_y)

        img2 = np.zeros([length, width, band], dtype=np.float32)
        for j in range(0, width):
            for i in range(0, length):
                a1 = length - size - i
                a2 = width - ysize - j
                if a1 < 0:
                    size2 = size + a1
                else:
                    size2 = size
                if a2 < 0:
                    size3 = size + a2
                else:
                    size3 = size
                size4 = min(size2, size3)
                blur_kernel = np.full(size, 1 / float(size4), dtype=np.float32)

                for k in range(0, size4):
                    ii = i + blur_kernel_x[k]
                    jj = j + blur_kernel_y[k]
                    img2[i, j, :] = img2[i, j, :] + \
                        img[ii, jj] * blur_kernel[k]
    return img2


def image_multi_process(func, width_n, height_n, out_band, out_dtype, inArray, args):

    ksize = inArray.shape[0]
    jsize = inArray.shape[1]
    isize = inArray.shape[2]
    ijksize = ksize, jsize, isize
    print(ijksize)

    split_data = split(inArray, width_n, height_n, ksize)
    jsize2 = split_data[0].shape[1]
    isize2 = split_data[0].shape[2]
    ijksize2 = out_band, jsize2, isize2
    ijksize3 = out_band, jsize, isize
    data_pre = []

    [data_pre.append(np.zeros(ijksize2, dtype=out_dtype))
     for n in range(width_n * height_n)]

    del inArray
    # --------calculation------------

    pipes = []
    jobs = []
    n3 = -1
    # octave number (octave number * core number = process number)
    for m in range(1, height_n + 1):

        for n in range(1, width_n + 1):  # core number
            n1 = width_n * (m - 1) + n - 1
            get_rev, send_rev = mp.Pipe(duplex=False)
            job = mp.Process(
                target=func,
                args=(
                    split_data[n1],
                    *args,
                    send_rev,
                ),
            )
            job.daemon = True
            jobs.append(job)
            pipes.append(get_rev)
            job.start()

            print(str(n1) + " " + str(n) + " " + str(m))

        for n in range(1, width_n + 1):
            n1 = width_n * (m - 1) + n - 1
            n3 = n3 + 1
            data_pre[n1] = pipes[n3].recv()

        [job.join() for job in jobs]

    data_o = reunion(data_pre, width_n, height_n, ijksize3)

    return data_o
