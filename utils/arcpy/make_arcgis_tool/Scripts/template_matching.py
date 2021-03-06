#!/usr/bin/env python
import sys
import numpy as np
import os
import respack
import arcpy
import scipy
import time
import ConfigParser
import subprocess
#-----------initialize------------
fort_path=os.path.split(os.path.split(respack.__file__)[0])[0]

tool_path=os.path.split(fort_path)[0]

arcpy.gp.overwriteOutput=True

#command = ["C:\\ProgramData\Anaconda3\\envs\\tecogan\\python.exe","D:\\deeplearning\\TecoGAN-master\\TecoGAN-master\\runGan.py","1"]
#subprocess.call(command,cwd="D:\\deeplearning\\TecoGAN-master\\TecoGAN-master")

inLayer = arcpy.GetParameterAsText(0)
outLayer1 = arcpy.GetParameterAsText(1)
outLayer2 = arcpy.GetParameterAsText(2)
outLayer3 = arcpy.GetParameterAsText(3)
outLayer4 = arcpy.GetParameterAsText(4)
tempDir = arcpy.GetParameterAsText(5)
tempImg = arcpy.GetParameterAsText(6)
mode1 = arcpy.GetParameter(7)
mode2 = arcpy.GetParameter(8)

bitNum = arcpy.GetParameter(9)
bandB = arcpy.GetParameter(10) -1
bandG = arcpy.GetParameter(11) - 1
bandR = arcpy.GetParameter(12) - 1
bandN = arcpy.GetParameter(13) - 1

bitNumt = arcpy.GetParameter(14)
bandBt = arcpy.GetParameter(15) - 1
bandGt = arcpy.GetParameter(16) - 1
bandRt = arcpy.GetParameter(17) - 1
bandNt = arcpy.GetParameter(18) - 1

cdepth=2**bitNum
cdeptht=2**bitNumt

arcpy.AddMessage("inLayer = {}".format(inLayer))
arcpy.AddMessage("outLayer1 = {}".format(outLayer1))
arcpy.AddMessage("outLayer2 = {}".format(outLayer2))
arcpy.AddMessage("outLayer3 = {}".format(outLayer3))
arcpy.AddMessage("outLayer4 = {}".format(outLayer4))
arcpy.AddMessage("tempDir = {}".format(tempDir))
arcpy.AddMessage("tempImg = {}".format(tempImg))
arcpy.AddMessage("mode1 = {}".format(mode1))
arcpy.AddMessage("mode2 = {}".format(mode2))


#--------------set parameters-------------------
if mode2:   #for high resolution
 parameters='HIGH_RESOLUTION'
else:       #for low resolution
 parameters='LOW_RESOLUTION' 
#------------------------------------------------
config=ConfigParser.ConfigParser()
inifile=os.path.join(tool_path,"Scripts\\template_matching.ini")
config.read(inifile)
ndvi_thres=config.get(parameters,'ndvi_thres')
confidence_thres=config.get(parameters,'confidence_thres')
scale_thres=config.get(parameters,'scale_thres')
bbs_thres=config.get(parameters,'bbs_thres')

arcpy.AddMessage("thres1 = {}".format(ndvi_thres))
arcpy.AddMessage("thres2 = {}".format(confidence_thres))
arcpy.AddMessage("thres3 = {}".format(scale_thres))
arcpy.AddMessage("thres4 = {}".format(bbs_thres))  
#---------------------------------------------


ltfile=os.path.join(fort_path,"src\\lookuptable.txt")

if mode2:
 width=17
 length=17
else:
 width=17
 length=17

if mode1:
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

# arrTempImgMod=np.zeros((1,ksize,length,width))
# arrTempImgMod=arrTempImgMod.astype(np.float32)
# arrTempImgRes=respack.template_matching_hr.resize(arrTempImg,ksize,length,width)
# arrTempImgMod[0,:,:,:]=arrTempImgRes[:,:,:]
# arrTempImgMod[:,0:4,:,:]=arrTempImgMod[:,0:4,:,:]/float(cdeptht) 
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
   
   arrTempImgRes=respack.template_matching_hr.resize(arrTempImg,ksize,length,width)    
   arrTempImgMod[num,:,:,:]=arrTempImgRes
   num=1+num
   
# arrTempImgMod=arrTempImg.astype(np.float32)
 arrTempImgMod[:,0:4,:,:]=arrTempImgMod[:,0:4,:,:]/float(cdeptht)
 
#---------------------------------------------
  
if mode2:
 template_matching=respack.template_matching_hr.landslide_detection
else:
 template_matching=respack.template_matching_lr.landslide_detection
 
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
data_o1=template_matching(inArray,arrTempImgMod,ltfile,ndvi_thres,bbs_thres,confidence_thres,scale_thres)
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
outRaster1 = arcpy.NumPyArrayToRaster(data_o1[0,:,:],llc,cx,cy,0)
arcpy.DefineProjection_management(outRaster1,sr)
outRaster2 = arcpy.NumPyArrayToRaster(data_o1[1,:,:],llc,cx,cy,0)
arcpy.DefineProjection_management(outRaster2,sr)
outRaster3 = arcpy.NumPyArrayToRaster(data_o1[2,:,:],llc,cx,cy,0)
arcpy.DefineProjection_management(outRaster3,sr)
outRaster4 = arcpy.NumPyArrayToRaster(data_o1[3,:,:],llc,cx,cy,0)
arcpy.DefineProjection_management(outRaster4,sr)
#-------test--------------
#outRaster2=arcpy.MakeRasterLayer_management(outRaster1, "outRaster2")
#arcpy.SetParameterAsText(1,outRaster2)
#------------------------
#arcpy.MakeRasterLayer_management(outRaster1, outLayer1)
#arcpy.MakeRasterLayer_management(outRaster2, outLayer2)

toutLayer1=arcpy.MakeRasterLayer_management(outRaster1,'tempLayer1')
toutLayer2=arcpy.MakeRasterLayer_management(outRaster2,'tempLayer2')

arcpy.CopyRaster_management(toutLayer1,outLayer1)
arcpy.CopyRaster_management(toutLayer2,outLayer2)
#outRaster2.save(toutLayer2)

arcpy.RasterToPolygon_conversion(outRaster3,outLayer3,"SIMPLIFY","VALUE")
arcpy.RasterToPolygon_conversion(outRaster4,outLayer4,"SIMPLIFY","VALUE")

#arcpy.SetParameterAsText(1,outLayer1)
#arcpy.SetParameterAsText(2,outLayer2)
#arcpy.SetParameterAsText(3,outLayer3)
#arcpy.SetParameterAsText(4,outLayer4)

endTime = float(time.time())
arcpy.AddMessage(endTime-startTime)
arcpy.AddMessage("Finish")
