#cythonize -i geometry.pyx, not cythonize -b. the latter results in memory error.
import numpy as np
import cython
cimport numpy as np
from cython import boundscheck, wraparound
from libc.math cimport sqrt,sin,cos,pi,acos

np.import_array()
DTYPE = np.float32
ctypedef np.float32_t DTYPE_t

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef map_2D(object func, np.ndarray[DTYPE_t, ndim=3] inArray):
    cdef int width = inArray.shape[1]
    cdef int height = inArray.shape[2]
    cdef int band = inArray.shape[0]
    cdef np.ndarray[DTYPE_t, ndim=2] outArray = np.empty([width, height], dtype=DTYPE)
    cdef int i,j
    for i in range(width):
        for j in range(height):
            outArray[i,j]=func(*[inArray[k,i,j] for k in range(band)])
    return outArray




