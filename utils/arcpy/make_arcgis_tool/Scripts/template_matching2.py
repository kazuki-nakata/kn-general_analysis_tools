#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import numpy as np
import os
import respack
import arcpy
import scipy
import time
import ConfigParser
import subprocess
from scipy.stats import truncnorm
import copy

def trunc_norm_distribution(parameters):
    a, b = (parameters[2] - parameters[0]) / parameters[1], (parameters[3] - parameters[0]) / parameters[1]
    x = np.arange(0, 10, 0.05)/10
    y = truncnorm.pdf(x, a, b, loc=parameters[0],scale=parameters[1])
    return y/sum(y)

def difference(p, q):
    kl = np.sum(abs(p - q)) *  1/2 
    return kl

def calc_gradient(num,parameters,parameter_size,used_func,base,interval):
    
    for i in range(0,parameter_size-2):
#        print(i)
        if num == i:
            f0=used_func(parameters=parameters,p=base)       
            parameters[i] = parameters[i]+ interval[i]
            f1=used_func(parameters=parameters,p=base)
            grad = (f1-f0)/interval[i]
            break    
    return grad

def calc_difference(parameters,p):
    q=trunc_norm_distribution(parameters)#normal_distribution(parameters)
    return difference(q,p)

def gradient_discent(used_pdf, parameters, base, parameter_size, learning_rate, max_epoch, eps):
    cond = eps + 10
    epoch = 0
    interval = np.zeros(parameter_size)
    interval[0:parameter_size] = 0.001
    used_func=calc_difference
    z0=used_func(parameters,base)
    epoch = 0
    while cond > eps and epoch < max_epoch:
        parameters_old = copy.copy(parameters)
        for i in range(0,parameter_size-2):
#            print("epoch",i,epoch)

            parameters[i] = parameters[i] - learning_rate * calc_gradient(i,parameters_old,parameter_size,used_func,base,interval)

        z1=used_func(parameters,base)
        cond = abs(z1-z0)
        z0=copy.copy(z1)
        
        epoch = epoch + 1

    parameters_new = parameters

    return parameters_new

#-----------initialize------------
fort_path=os.path.split(os.path.split(respack.__file__)[0])[0]
tool_path=os.path.split(fort_path)[0]
arcpy.gp.overwriteOutput=True

inLayer = arcpy.GetParameterAsText(0)

outLayer1 = arcpy.GetParameterAsText(1)
tempDir = arcpy.GetParameterAsText(2)
tempImg = arcpy.GetParameterAsText(3)
temp_mode = arcpy.GetParameter(4)
para_mode = arcpy.GetParameter(5)

sigma_coff = arcpy.GetParameter(24)

ndvi_thres=arcpy.GetParameter(7)
confidence_thres=arcpy.GetParameter(8)
scale_thres=arcpy.GetParameter(9)
sim_bbs_thres=arcpy.GetParameter(10)
sim_ndvi_thres=arcpy.GetParameter(11)
sim_br_thres=arcpy.GetParameter(12)
sim_gsi_thres=arcpy.GetParameter(13)

outLayer2 = arcpy.GetParameterAsText(6)

bitNum = arcpy.GetParameter(14)
bandB = arcpy.GetParameter(15) -1
bandG = arcpy.GetParameter(16) - 1
bandR = arcpy.GetParameter(17) - 1
bandN = arcpy.GetParameter(18) - 1

bitNumt = arcpy.GetParameter(19)
bandBt = arcpy.GetParameter(20) - 1
bandGt = arcpy.GetParameter(21) - 1
bandRt = arcpy.GetParameter(22) - 1
bandNt = arcpy.GetParameter(23) - 1

cdepth=2**bitNum
cdeptht=2**bitNumt

arcpy.AddMessage(bitNum)
arcpy.AddMessage(bandB)
arcpy.AddMessage(bandG)
arcpy.AddMessage(bandR)
arcpy.AddMessage(bandN)
arcpy.AddMessage(inLayer)
arcpy.AddMessage(inLayer)
arcpy.AddMessage(outLayer1)
arcpy.AddMessage(outLayer2)
arcpy.AddMessage(tempDir)
arcpy.AddMessage(tempImg)
arcpy.AddMessage(cdepth)
arcpy.AddMessage("temp_mode = {}".format(temp_mode))
arcpy.AddMessage("para_mode = {}".format(para_mode))
arcpy.AddMessage("ndvi_thres = {}".format(ndvi_thres))
arcpy.AddMessage("conf_thres = {}".format(confidence_thres))
arcpy.AddMessage("scale_thres = {}".format(scale_thres))
arcpy.AddMessage("sim_bbs_thres = {}".format(sim_bbs_thres))
arcpy.AddMessage("sim_br_thres = {}".format(sim_br_thres))

if temp_mode:
    arcpy.AddMessage("temp_mode:True")
else:
    arcpy.AddMessage("temp_mode:False")
    
if para_mode:
    arcpy.AddMessage("para_mode:True")
else:
    arcpy.AddMessage("para_mode:False")   

if len(outLayer2) != 0:
    arcpy.AddMessage(outLayer2)
else:
    arcpy.AddMessage("outLayer2:None")

if sim_ndvi_thres is not None:
    arcpy.AddMessage("sim_ndvi_thres = {}".format(sim_ndvi_thres))
else:
    arcpy.AddMessage("sim_ndvi_thres:None")
    sim_ndvi_thres = -1

if sim_br_thres is not None:
    arcpy.AddMessage("sim_brightness_thres = {}".format(sim_br_thres))
else:
    arcpy.AddMessage("sim_brightness_thres:None")
    sim_br_thres = -1    

if sim_gsi_thres is not None:
    arcpy.AddMessage("sim_gsi_thres = {}".format(sim_gsi_thres))
else:
    arcpy.AddMessage("sim_gsi_thres:None")
    sim_gsi_thres = -1  
    
width=17
length=17

if temp_mode:   
    rasTempImg = arcpy.Raster(tempImg)
    arrTempImg = arcpy.RasterToNumPyArray(rasTempImg)
    ksize=arrTempImg.shape[0]
    jsize=arrTempImg.shape[1]
    isize=arrTempImg.shape[2]

    if bandBt >= 0:
        arrTempImg[[bandBt],:,:],arrTempImg[[0],:,:]=arrTempImg[[0],:,:],arrTempImg[[bandBt],:,:]

    if bandGt >= 1:
        arrTempImg[[bandGt],:,:],arrTempImg[[1],:,:]=arrTempImg[[1],:,:],arrTempImg[[bandGt],:,:]

    if bandRt >= 2:
        arrTempImg[[bandRt],:,:],arrTempImg[[2],:,:]=arrTempImg[[2],:,:],arrTempImg[[bandRt],:,:]

    if bandNt >= 3:
        arrTempImg[[bandNt],:,:],arrTempImg[[3],:,:]=arrTempImg[[3],:,:],arrTempImg[[bandNt],:,:]

    arrTempImgMod=np.zeros((1,ksize,jsize,isize))
    arrTempImgMod=arrTempImgMod.astype(np.float32)
    arrTempImgMod[0,:,:,:]=arrTempImg[:,:,:]
    arrTempImgMod[:,0:4,:,:]=arrTempImgMod[:,0:4,:,:]#/float(cdeptht)

else:

    file_list = os.listdir(tempDir)

    num2=0
    for file_name in file_list:
        root, ext = os.path.splitext(file_name)
        if ext == u'.tif':
            num2=num2+1

    num=0
    for file_name in file_list:
        root, ext = os.path.splitext(file_name)
  
        if ext == u'.tif':
            tempImg = tempDir + '\\' + file_name
            arcpy.AddMessage("InputTemplate = {}".format(tempImg))
            rasTempImg = arcpy.Raster(tempImg)
            arrTempImg = arcpy.RasterToNumPyArray(rasTempImg)
            if bandBt >= 0:
                arrTempImg[[bandBt],:,:],arrTempImg[[0],:,:]=arrTempImg[[0],:,:],arrTempImg[[bandBt],:,:]

            if bandGt >= 1:
                arrTempImg[[bandGt],:,:],arrTempImg[[1],:,:]=arrTempImg[[1],:,:],arrTempImg[[bandGt],:,:]

            if bandRt >= 2:
                arrTempImg[[bandRt],:,:],arrTempImg[[2],:,:]=arrTempImg[[2],:,:],arrTempImg[[bandRt],:,:]

            if bandNt >= 3:
                arrTempImg[[bandNt],:,:],arrTempImg[[3],:,:]=arrTempImg[[3],:,:],arrTempImg[[bandNt],:,:]

            if num ==0:
                ksize=arrTempImg.shape[0]
                jsize=arrTempImg.shape[1]
                isize=arrTempImg.shape[2]

                arrTempImgMod=np.zeros((num2,ksize,length,width))
                arrTempImgMod=arrTempImgMod.astype(np.float32)

                print str(width) + '-' + str(length)
   
            arrTempImgRes=respack.advanced_process.resize(arrTempImg,ksize,length,width)    
            arrTempImgMod[num,:,:,:]=arrTempImgRes
            num=1+num
   
    arrTempImgMod[:,0:4,:,:]=arrTempImgMod[:,0:4,:,:]#/float(cdeptht)

desc = arcpy.Describe(inLayer)
layerPath = desc.catalogPath
descPath = desc.path
outputPath = os.path.join(descPath, "output.shp")

#----------------------------------------
arcpy.AddMessage("Convert Raster into Array")
startTime = float(time.time())
inRaster = arcpy.Raster(inLayer)
inArray = arcpy.RasterToNumPyArray(inRaster)

if bandB >= 0:
    inArray[[bandB],:,:],inArray[[0],:,:]=inArray[[0],:,:],inArray[[bandB],:,:]

if bandG >= 1:
    inArray[[bandG],:,:],inArray[[1],:,:]=inArray[[1],:,:],inArray[[bandG],:,:]

if bandR >= 2:
    inArray[[bandR],:,:],inArray[[2],:,:]=inArray[[2],:,:],inArray[[bandR],:,:]

if bandN >= 3:
    inArray[[bandN],:,:],inArray[[3],:,:]=inArray[[3],:,:],inArray[[bandN],:,:]

inArray=inArray.astype(np.float32)
inArray[0:4,:,:]=inArray[0:4,:,:]/float(cdepth)

ksize=inArray.shape[0]
jsize=inArray.shape[1]
isize=inArray.shape[2]

ijksize2=4,jsize,isize
ijksize=ksize,jsize,isize
ijsize=jsize,isize
arcpy.AddMessage(ijksize)
data_o1=np.zeros(ijksize2,dtype=np.int32)
endTime = float(time.time())
arcpy.AddMessage(endTime-startTime)

#--------calculation------------
arcpy.AddMessage("Calculation")
arcpy.AddMessage(arrTempImgMod.shape)
startTime = float(time.time())

#------test----------------------
cx= inRaster.meanCellWidth
cy= inRaster.meanCellHeight
mx = inRaster.extent.XMin
my = inRaster.extent.YMin
llc=arcpy.Point(mx,my)
sr=inRaster.spatialReference

if para_mode:
    hist_par = 0.99
    bin_width=0.005        
    bin_number = round(1/bin_width)
#    arcpy.AddMessage(bin_number)     

    proc_histgram = respack.advanced_process.calc_histgram_parameter
    hmax,hmode,pdf = proc_histgram(inArray,bin_width,bin_number,hist_par)
    arcpy.AddMessage(hmode)
    arcpy.AddMessage(hmax)  
    
    amin=hmode - (hmax-hmode)
    amax=hmode + (hmax-hmode)

    _,hsize=pdf.shape

    pdf2 = np.zeros((2,hsize),dtype=np.float)
    pdf2[0,:] = pdf[0,:]

    for i in range(0,hsize):
        if pdf[0,i] < amin or pdf[0,i] > amax:
            pdf2[1,i]=0
        else:
            pdf2[1,i]=pdf[1,i]

    a=sum(pdf2[1,]) 
    pdf2[1,]=pdf2[1,]/a

    x = np.arange(0, 10, 0.05)/10

    #-------------------------------------------------
    used_pdf = trunc_norm_distribution
    learning_rate = 0.0001
    max_epoch = 10000
    parameter_size = 4
    eps = 0.000001

    parameters = np.zeros(parameter_size)
    parameters[:] = hmode,0.1,amin,amax

    parameters_new = gradient_discent(used_pdf, parameters, pdf2[1], parameter_size, learning_rate, max_epoch, eps)
    ndvi_thres = parameters_new[0] - sigma_coff * parameters_new[1]
    arcpy.AddMessage(parameters_new)
    arcpy.AddMessage(ndvi_thres)
    
proc_mask = respack.advanced_process.make_maskdata2
mask = proc_mask(inArray,ndvi_thres)
arcpy.AddMessage("Finish mask process")
if len(outLayer2) != 0:
    outRaster2 = arcpy.NumPyArrayToRaster(mask[0,:,:],llc,cx,cy,999)
    arcpy.DefineProjection_management(outRaster2,sr)
    toutLayer2=arcpy.MakeRasterLayer_management(outRaster2,'tempLayer1')
    arcpy.CopyRaster_management(toutLayer2,outLayer2)

proc_scale = respack.advanced_process.scale_selection
scale = proc_scale(inArray, confidence_thres, scale_thres)
arcpy.AddMessage("Finish scale process")
# # outRaster2 = arcpy.NumPyArrayToRaster(scale[1,:,:],llc,cx,cy,999)
# # arcpy.DefineProjection_management(outRaster2,sr)
# # toutLayer2=arcpy.MakeRasterLayer_management(outRaster2,'tempLayer1')
# # arcpy.CopyRaster_management(toutLayer2,outLayer2)

# # outRaster2 = arcpy.NumPyArrayToRaster(scale[0,:,:],llc,cx,cy,999)
# # arcpy.DefineProjection_management(outRaster2,sr)
# # toutLayer2=arcpy.MakeRasterLayer_management(outRaster2,'tempLayer1')
# # arcpy.CopyRaster_management(toutLayer2,outLayer3)

proc_window = respack.advanced_process.create_window_candidate
window = proc_window(scale[1,:,:],scale[0,:,:],mask[1,:,:], confidence_thres,scale_thres)
arcpy.AddMessage("Finish window process")

proc_orient = respack.advanced_process.orientation_selection
orient = proc_orient(inArray,scale[0,:,:],window)
arcpy.AddMessage("Finish orient process")

scan = 3 

proc_match = respack.advanced_process.template_matching2
data_o1 = proc_match(inArray[0:4,:,:],scale[0,:,:],scale[1,:,:],orient,window,arrTempImgMod,sim_bbs_thres,sim_br_thres,sim_ndvi_thres,sim_gsi_thres,scan)
arcpy.AddMessage("Finish matching process")
data_o2 = data_o1[1,:,:]
data_o2 = data_o2.astype(np.int32)
outRaster2 = arcpy.NumPyArrayToRaster(data_o2,llc,cx,cy,0)
arcpy.DefineProjection_management(outRaster2,sr)
arcpy.RasterToPolygon_conversion(outRaster2,outLayer1,"SIMPLIFY","VALUE")

endTime = float(time.time())
arcpy.AddMessage(endTime-startTime)


#-------output---------
# arcpy.AddMessage("Make the result layers")
# startTime = float(time.time())
# cx= inRaster.meanCellWidth
# cy= inRaster.meanCellHeight
# mx = inRaster.extent.XMin
# my = inRaster.extent.YMin

# llc=arcpy.Point(mx,my)
# sr=inRaster.spatialReference
# outRaster1 = arcpy.NumPyArrayToRaster(data_o1[0,:,:],llc,cx,cy,0)
# arcpy.DefineProjection_management(outRaster1,sr)
# outRaster2 = arcpy.NumPyArrayToRaster(data_o1[1,:,:],llc,cx,cy,0)
# arcpy.DefineProjection_management(outRaster2,sr)
# outRaster3 = arcpy.NumPyArrayToRaster(data_o1[2,:,:],llc,cx,cy,0)
# arcpy.DefineProjection_management(outRaster3,sr)
# outRaster4 = arcpy.NumPyArrayToRaster(data_o1[3,:,:],llc,cx,cy,0)
# arcpy.DefineProjection_management(outRaster4,sr)

# toutLayer1=arcpy.MakeRasterLayer_management(outRaster1,'tempLayer1')
# toutLayer2=arcpy.MakeRasterLayer_management(outRaster2,'tempLayer2')

# arcpy.CopyRaster_management(toutLayer1,outLayer1)
# arcpy.CopyRaster_management(toutLayer2,outLayer2)

# arcpy.RasterToPolygon_conversion(outRaster3,outLayer3,"SIMPLIFY","VALUE")
# arcpy.RasterToPolygon_conversion(outRaster4,outLayer4,"SIMPLIFY","VALUE")

# endTime = float(time.time())
# arcpy.AddMessage(endTime-startTime)
arcpy.AddMessage("Finish")
