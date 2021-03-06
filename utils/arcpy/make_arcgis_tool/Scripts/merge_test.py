#!/usr/bin/env python
import sys
import numpy as np
import os
import respack
import arcpy
import time
#-----------initialize------------
fort_path=os.path.split(os.path.split(respack.__file__)[0])[0]
ltfile=os.path.join(fort_path,"src\\lookuptable.txt")
arcpy.gp.overwriteOutput=True
#inLayer = "D:\\JR\\Composite_IMG_SPOT6_PMS_201807220125106_ORT_SPOT6_20180723_07162119f5m1gt5po9n_1_R1C11.tif"
#inLayer = "D:\\JR\\Composite_IMG_SPOT6_PMS_201807220125106_ORT_SPOT6_20180723_0717051k90t929yt1gx_1_R1C11.tif"
inLayer = "D:\\RESTEC_ARCGIS_TOOLS\\Demo\\template_matching_worldview3.tif"
#inLayer = "D:\\JR\\test\\yaki.tif"
#inLayer = "D:\\RESTEC_ARCGIS_TOOLS\\Demo\\template_matching_landsat8.tif"
outLayer ="test.tif"
tempDir = r"D:\RESTEC_ARCGIS_TOOLS\Fort_library\share\LC_N36E137.tif"
#tempImg = "D:\\RESTEC_ARCGIS_TOOLS\\Fort_library\\share\\template_image\\low_resolution\\template.tif"
#tempDir = "D:\\JR\\template"
#tempImg = "D:\\JR\\testtemp.tif"
mode1 = False
mode2 = True

if(mode2):
 ndvi_thres=0.4  #0.2-0.4
 bbs_thres=0.3   #0.2-0.4
else:
 ndvi_thres=0.3  #0.2-0.4
 bbs_thres=0.4   #0.2-0.4


bitNum = 11
bandB = 1 - 1
bandG = 2 - 1
bandR = 3 - 1
bandN = 4 - 1

bitNumt = 11
bandBt = 1 - 1
bandGt = 2 - 1
bandRt = 3 - 1
bandNt = 4 - 1

cdepth=2**bitNum
cdeptht=2**bitNumt

if(mode2):
 width=17
 length=17
else:
 width=17
 length=17

if(mode1):
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
 print arrTempImgMod[0,:,9,9] 
 arrTempImgMod[:,0:4,:,:]=arrTempImgMod[:,0:4,:,:]/float(cdeptht)
 print arrTempImgMod[0,:,9,9]
# print arraTempImgMod
else:

 file_list = os.listdir(tempDir)
 
 num2=0
 for file_name in file_list:
  root, ext = os.path.splitext(file_name)
  if ext == u'.tif':
   num2=num2+1

 print 'Number of template image =' + str(num2)
 num=0   
 for file_name in file_list:
  root, ext = os.path.splitext(file_name)
  print root + ext
  if ext == u'.tif':
   tempImg = tempDir + '\\' + file_name
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
    print arrTempImgMod.shape[0]
    
   print str(width) + '-' + str(length)
   
   arrTempImgRes=respack.template_matching_hr.resize(arrTempImg,ksize,length,width)
   
   arrTempImgMod[num,:,:,:]=arrTempImgRes
   num=1+num
   
 #arrTempImgMod=arrTempImgMod.astype(np.float32) 
 print arrTempImgMod[0,:,9,9]
 arrTempImgMod[:,0:4,:,:]=arrTempImgMod[:,0:4,:,:]/float(cdeptht)
 print arrTempImgMod[0,:,9,9]
 
#---------------------------------------------
  
if(mode2):
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
#print inArray[:,500,500]
inArray=inArray.astype(np.float32)
inArray[0:4,:,:]=inArray[0:4,:,:]/float(cdepth)

#print inArray[:,500,500]

ksize=inArray.shape[0]
jsize=inArray.shape[1]
isize=inArray.shape[2]

ijksize=ksize,jsize,isize
ijsize=jsize,isize
arcpy.AddMessage(ijksize)
data_o1=np.zeros(ijsize,dtype=np.int32)
endTime = float(time.time())
arcpy.AddMessage(endTime-startTime)

print 'finish'

#--------calculation------------
arcpy.AddMessage("Calculation")
startTime = float(time.time())
data_o1=template_matching(inArray,arrTempImgMod,ltfile,ndvi_thres,bbs_thres)
endTime = float(time.time())
arcpy.AddMessage(endTime-startTime)

#-------output---------
#arcpy.AddMessage("Make the result layers")
#startTime = float(time.time())
#cx= inRaster.meanCellWidth
#cy= inRaster.meanCellHeight
#mx = inRaster.extent.XMin
#my = inRaster.extent.YMin

#llc=arcpy.Point(mx,my)
#sr=inRaster.spatialReference
#outRaster1 = arcpy.NumPyArrayToRaster(data_o1,llc,cx,cy,0)
#arcpy.DefineProjection_management(outRaster1,sr)
#outRaster1.save("D:\\RESTEC_ARCGIS_TOOLS\\Demo\\test5.tif")
#print cx,cy,mx,my


#-------test--------------
#outRaster2=arcpy.MakeRasterLayer_management(outRaster1, "test.tif")
#arcpy.SetParameterAsText(1,outRaster2)
#------------------------

#arcpy.RasterToPolygon_conversion(outRaster1,'D:\\RESTEC_ARCGIS_TOOLS\\Demo\\test.shp',"SIMPLIFY","VALUE")
#arcpy.SetParameterAsText(1,outLayer)
#endTime = float(time.time())
#arcpy.AddMessage(endTime-startTime)
#arcpy.AddMessage("Finish")
