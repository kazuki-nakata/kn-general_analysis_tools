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
#    print parameters
#    print y/sum(y)
    return y/sum(y)

def difference(p, q):
#    print p-q
    kl = np.sum(abs(p - q)) * 1/2 
#    print kl
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
#    print q
    return difference(q,p)

def gradient_discent(used_pdf, parameters, base, parameter_size, learning_rate, max_epoch, eps):
    cond = eps + 10
    epoch = 0
    interval = np.zeros(parameter_size)
    interval[0:parameter_size] = 0.001
    used_func=calc_difference
#    print parameters,base
    z0=used_func(parameters,base)
#    print z0
#    print "hi"
    epoch = 0
    while cond > eps and epoch < max_epoch:
        parameters_old = copy.copy(parameters)
        for i in range(0,parameter_size-2):
#            print("epoch",i,epoch)

            parameters[i] = parameters[i] - learning_rate * calc_gradient(i,parameters_old,parameter_size,used_func,base,interval)

        z1=used_func(parameters,base)
        cond = abs(z1-z0)
        z0=copy.copy(z1)
#        print z1,z0
        epoch = epoch + 1
#        print parameters, cond
    parameters_new = parameters

    return parameters_new


fort_path=os.path.split(os.path.split(respack.__file__)[0])[0]
ltfile=os.path.join(fort_path,"src\\lookuptable.txt")
arcpy.gp.overwriteOutput=True


inLayer = arcpy.GetParameterAsText(0)

outLayer1 = arcpy.GetParameterAsText(1)
tempDir = arcpy.GetParameterAsText(2)
tempImg = arcpy.GetParameterAsText(3)
temp_mode = arcpy.GetParameter(4)
para_mode = arcpy.GetParameter(5)

# ndvi_thres=arcpy.GetParameter(7)
# confidence_thres=arcpy.GetParameter(8)
# scale_thres=arcpy.GetParameter(9)
# sim_bbs_thres=arcpy.GetParameter(10)
# sim_ndvi_thres=arcpy.GetParameter(11)
# sim_brightness_thres=arcpy.GetParameter(12)
# sim_gsi_thres=arcpy.GetParameter(13)

# outLayer2 = arcpy.GetParameterAsText(6)

# bitNum = arcpy.GetParameter(14)
# bandB = arcpy.GetParameter(15) -1
# bandG = arcpy.GetParameter(16) - 1
# bandR = arcpy.GetParameter(17) - 1
# bandN = arcpy.GetParameter(18) - 1

# bitNumt = arcpy.GetParameter(19)
# bandBt = arcpy.GetParameter(20) - 1
# bandGt = arcpy.GetParameter(21) - 1
# bandRt = arcpy.GetParameter(22) - 1
# bandNt = arcpy.GetParameter(23) - 1
# inLayer = r"D:\JR\2019\ツールテスト\マスク\result_shape.tif"
# outLayer1 =r"D:\JR\2019\ツールテスト\抽出\test.shp"
# outLayer2 = r"D:\JR\2019\ツールテスト\抽出\test.tif"
# outLayer3 = r"D:\JR\2019\ツールテスト\抽出\test2.tif"
# tempDir = r"D:\RESTEC_ARCGIS_TOOLS_DEV\Fort_library\share\template_image"

temp_mode = False
para_mode = False

ndvi_thres=0.55
confidence_thres=0.01
scale_thres=3
sim_bbs_thres=0.25
sim_ndvi_thres=0
sim_brightness_thres=0
sim_gsi_thres=0
scan = 3 


bitNum = 12
bandB = 3 - 1
bandG = 2 - 1 
bandR = 1 - 1
bandN = 4 - 1

bitNumt = 11
bandBt = 1 - 1
bandGt = 2 - 1
bandRt = 3 - 1
bandNt = 4 - 1


cdepth=2**bitNum
cdeptht=2**bitNumt


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
 arrTempImgMod[:,0:4,:,:]=arrTempImgMod[:,0:4,:,:]/float(cdeptht)

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
   

 arrTempImgMod[:,0:4,:,:]=arrTempImgMod[:,0:4,:,:]/float(cdeptht)

desc = arcpy.Describe(inLayer)
layerPath = desc.catalogPath
descPath = desc.path
outputPath = os.path.join(descPath, "output.shp")

#----------------------------------------

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


para_mode = True
if para_mode:
    hist_par = 0.99
    bin_width=0.005        
    bin_number = round(1/bin_width)


    proc_histgram = respack.advanced_process.calc_histgram_parameter
    hmax,hmode,pdf = proc_histgram(inArray,bin_width,bin_number,hist_par)

    amin=hmode - (hmax-hmode)
    amax=hmode + (hmax-hmode)
    arcpy.AddMessage(amax)
    _,hsize=pdf.shape
    #    print pdf
    pdf2 = np.zeros((2,hsize),dtype=np.float)
    pdf2[0,:] = pdf[0,:]

    for i in range(0,hsize):
        if pdf[0,i] < amin or pdf[0,i] > amax:
            pdf2[1,i]=0
        else:
            pdf2[1,i]=pdf[1,i]

    a=sum(pdf2[1,:]) 
    pdf2[1,:]=pdf2[1,:]/a

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
    ndvi_thres = parameters_new[0] - 1.5 * parameters_new[1]
    #    print parameters_new,hmode,hmax,ndvi_thres
    arcpy.AddMessage(parameters_new)
    arcpy.AddMessage(ndvi_thres)    
# proc_mask = respack.advanced_process.make_maskdata2
# mask = proc_mask(inArray,ndvi_thres)
# arcpy.AddMessage("Finish mask process")
# outRaster2 = arcpy.NumPyArrayToRaster(mask[1,:,:],llc,cx,cy,999)
# arcpy.DefineProjection_management(outRaster2,sr)
# toutLayer2=arcpy.MakeRasterLayer_management(outRaster2,'tempLayer1')
# arcpy.CopyRaster_management(toutLayer2,outLayer2)

# outRaster2 = arcpy.NumPyArrayToRaster(mask[0,:,:],llc,cx,cy,999)
# arcpy.DefineProjection_management(outRaster2,sr)
# toutLayer2=arcpy.MakeRasterLayer_management(outRaster2,'tempLayer1')
# arcpy.CopyRaster_management(toutLayer2,outLayer3)

# proc_scale = respack.advanced_process.scale_selection
# scale = proc_scale(inArray, confidence_thres, scale_thres)
# arcpy.AddMessage("Finish scale process")
# outRaster2 = arcpy.NumPyArrayToRaster(scale[1,:,:],llc,cx,cy,999)
# arcpy.DefineProjection_management(outRaster2,sr)
# toutLayer2=arcpy.MakeRasterLayer_management(outRaster2,'tempLayer1')
# arcpy.CopyRaster_management(toutLayer2,outLayer2)

# outRaster2 = arcpy.NumPyArrayToRaster(scale[0,:,:],llc,cx,cy,999)
# arcpy.DefineProjection_management(outRaster2,sr)
# toutLayer2=arcpy.MakeRasterLayer_management(outRaster2,'tempLayer1')
# arcpy.CopyRaster_management(toutLayer2,outLayer3)

# proc_window = respack.advanced_process.create_window_candidate
# window = proc_window(scale[1,:,:],scale[0,:,:],mask[1,:,:], confidence_thres,scale_thres)
# arcpy.AddMessage("Finish window process")

# proc_orient = respack.advanced_process.orientation_selection
# orient = proc_orient(inArray,scale[0,:,:],window)
# arcpy.AddMessage("Finish orient process")

# proc_match = respack.advanced_process.template_matching2
# data_o1 = proc_match(inArray[0:4,:,:],scale[0,:,:],scale[1,:,:],orient,window,arrTempImgMod,sim_bbs_thres,scan)
# arcpy.AddMessage("Finish matching process")
# outRaster2 = arcpy.NumPyArrayToRaster(data_o1[1,:,:],llc,cx,cy,0)
# arcpy.DefineProjection_management(outRaster2,sr)

#toutLayer2=arcpy.MakeRasterLayer_management(outRaster2,'tempLayer1')
#arcpy.CopyRaster_management(toutLayer2,outLayer2)

# outRaster2 = arcpy.Raster(outLayer2)
# inArray = arcpy.RasterToNumPyArray(outRaster2)
# inArray = inArray.astype(np.int32)
# outRaster2 = arcpy.NumPyArrayToRaster(inArray,llc,cx,cy,0)
# arcpy.DefineProjection_management(outRaster2,sr)
# arcpy.RasterToPolygon_conversion(outRaster2,outLayer1,"SIMPLIFY","VALUE")

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
print "finish"