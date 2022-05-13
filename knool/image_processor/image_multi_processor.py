import multiprocessing as mp
from multiprocessing.connection import wait
import numpy as np
from general_image_treatment import image_split,image_reunion

class split_to_reunion():
    def main(self,func,width_n,height_n,out_band,out_dtype,inArray,args):

        ksize = inArray.shape[0]
        jsize = inArray.shape[1]
        isize = inArray.shape[2]
        ijksize = ksize, jsize, isize
        print(ijksize)
        ijsize = jsize, isize

        split_data = image_split(inArray, width_n, height_n, ksize)
        jsize2 = split_data[0].shape[1]
        isize2 = split_data[0].shape[2]
        ijksize2 = out_band,jsize2,isize2
        ijksize3 = out_band,jsize,isize
        data_pre = []

        [data_pre.append(np.zeros(ijksize2, dtype=out_dtype)) for n in range(width_n * height_n)]

        del inArray
        # --------calculation------------

        pipes = []
        jobs = []
        n3 = -1
        for m in range(1, height_n + 1):  # octave number (octave number * core number = process number)

            for n in range(1, width_n + 1):  # core number
                n1 = width_n * (m - 1) + n - 1
                get_rev, send_rev = mp.Pipe(duplex=False)
                job = mp.Process(target=func, args=(split_data[n1], *args, send_rev,))
                job.daemon = True
                jobs.append(job)
                pipes.append(get_rev)
                job.start()

                print(str(n1) + ' ' + str(n) + " " + str(m))

            for n in range(1, width_n + 1):
                n1 = width_n * (m - 1) + n - 1
                n3 = n3 + 1
                data_pre[n1] = pipes[n3].recv()

            [job.join() for job in jobs]

        data_o = image_reunion(data_pre, width_n, height_n, ijksize3)

        return data_o