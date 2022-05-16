import cython

cpdef calc_IST2(float bt1, float bt2, float senz):
    cdef float val
    if bt1<-240:
        val=bt1 + (bt1-bt2) + (bt1-bt2)
    elif bt1<-120:
        val=bt1 + (bt1-bt2) + (bt1-bt2)
    else:
        val=bt1 + (bt1-bt2) + (bt1-bt2)
    return val