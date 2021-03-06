#!/usr/bin/env python
import sys
import numpy as np
import os
import respack
import arcpy
import time
import multiprocessing as mp

def imageSplit(img,n):
    print img.shape[1]
    print img.shape[2]    
    pv_size = img.shape[1] // n
    ph_size = img.shape[2] // n
    v_size = img.shape[1] // pv_size * pv_size
    h_size = img.shape[2] // ph_size * ph_size
    img=img[0:2,:v_size, :h_size]
    out_img = []
    [out_img.extend(np.split(h_img, n, axis=2))
      for h_img in np.split(img, n, axis=1)]
    return out_img

def imageReunion(img,n,ijsize):
    out_img= np.zeros(ijsize,dtype=np.int32)
    h=img[0].shape[1]
    v=img[0].shape[0]    
    print h
    print v
    for vnum in range(n):
        for hnum in range(n):        
            n1 = vnum * n + hnum
            hm = h * hnum
            hl = h * (hnum + 1)
            vm = v * vnum
            vl = v * (vnum + 1)
            out_img[vm:vl,hm:hl]=img[n1][:,:]
    return out_img
    
def proc(img,send_rev):
    color_thres=1.2
    sigma_thres=80
    comp_thres=1 #0.35
    asp_thres=1 #0.8
    size_min_thres=60 #pixel    
    segmentation=respack.segmentation_rb.landslide_detection    
    print('process id:', os.getpid())
#    output=img[0,:,:]
    output=segmentation(img,sigma_thres,color_thres,comp_thres,asp_thres,size_min_thres)
    send_rev.send(output)

if __name__ == '__main__':

#-----------initialize------------
    #inLayer = "D:\\RESTEC_ARCGIS_TOOLS\\Demo\\segmentation_palsar2.tif"
    #segmentation=respack.segmentationRGB.landslide_detection
    inLayer = "D:\\segmentation\\FY30\\subset.tif"
    outLayer ="test.tif"
    arcpy.AddMessage("Convert Raster into Array")
    startTime = float(time.time())
    inRaster = arcpy.Raster(inLayer)
    inArray = arcpy.RasterToNumPyArray(inRaster)
    ksize=inArray.shape[0]
    jsize=inArray.shape[1]
    isize=inArray.shape[2]
    ijksize=ksize,jsize,isize
    ijksize2=2,jsize,isize
    ijsize=jsize,isize
    arcpy.AddMessage(ijksize)
    data_o= np.zeros(ijsize,dtype=np.int32)
    endTime = float(time.time())
    arcpy.AddMessage(endTime-startTime)

    split_data=imageSplit(inArray,4)
    jsize2=split_data[0].shape[1]
    isize2=split_data[0].shape[2]
    ijsize2=jsize2,isize2
    data_pre=[]

    [data_pre.append(np.zeros(ijsize2,dtype=np.int32)) for n in range(16)]
    
    del inArray
    #--------calculation------------
    arcpy.AddMessage("Calculation")
    startTime = float(time.time())

    pipes = []
    jobs = []
    n2 = -1
    
    for m in range(1,5): #octave number (octave number * core number = 16)
        for n in range(1,5): #core number
            n1 = 4 * (m - 1) + n - 1
            n2 = n2 + 1
            get_rev,send_rev = mp.Pipe(False)
            job = mp.Process(target=proc, args=(split_data[n1],send_rev))       
            jobs.append(job)
            pipes.append(get_rev)
            job.start()
            print str(n1) + ' ' + str(n2) + ' ' + str(n) + " " + str(m)
            data_pre[n1] = pipes[n2].recv()
        [data_pre[n1] for job in jobs]

    data_o=imageReunion(data_pre,4,ijsize)
    print str(data_o.shape[0]) + ',' + str(data_o.shape[1]) 
    endTime = float(time.time())
    arcpy.AddMessage(endTime-startTime)

    #-------output---------
    arcpy.AddMessage("Make the result layers")
    startTime = float(time.time())
    cx= inRaster.meanCellWidth
    cy= inRaster.meanCellHeight
    mx = inRaster.extent.XMin
    my = inRaster.extent.YMin

    llc=arcpy.Point(mx,my)
    sr=inRaster.spatialReference
    outRaster1 = arcpy.NumPyArrayToRaster(data_o,llc,cx,cy,0)
    arcpy.DefineProjection_management(outRaster1,sr)
    outRaster1.save("D:\\segmentation\\FY30\\reunion.tif")

    arcpy.RasterToPolygon_conversion(outRaster1,'D:\\segmentation\\FY30\\test.shp',"SIMPLIFY","VALUE")
    arcpy.SetParameterAsText(1,outLayer)
    endTime = float(time.time())
    arcpy.AddMessage(endTime-startTime)
    arcpy.AddMessage("Finish")
