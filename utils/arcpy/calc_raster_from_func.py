#This program is run only on the arcpy environment!

import sys
import numpy as np
import os
import arcpy
import time

#-----------initialize------------
arcpy.gp.overwriteOutput=True

#Only arcpy toolbox
# inLayer = arcpy.GetParameterAsText(0)
# outLayer = arcpy.GetParameterAsText(1)
# arcpy.AddMessage("inLayer = {}".format(inLayer))
# desc = arcpy.Describe(inLayer)
# layerPath = desc.catalogPath
# descPath = desc.path
# outputPath = os.path.join(descPath, "output.shp")

inLayer=r'inpath'
outLayer=r'outpath'
func='func_object'
#----------------------------------------
arcpy.AddMessage("Convert Raster into Array")
startTime = time.time()
inRaster = arcpy.Raster(inLayer)
inArray = arcpy.RasterToNumPyArray(inRaster)
ksize=inArray.shape[0]
jsize=inArray.shape[1]
isize=inArray.shape[2]
ijksize=ksize,jsize,isize
ijsize=jsize,isize
arcpy.AddMessage(ijksize)
data_o1=np.zeros(ijsize,dtype=np.int32)
endTime = float(time.time())
arcpy.AddMessage(endTime-startTime)

threshold = 0.6
#--------calculation------------
arcpy.AddMessage("Calculation")
startTime = float(time.time())
data_o1=func(inArray)
endTime = float(time.time())
arcpy.AddMessage(endTime-startTime)

#-------output---------
arcpy.AddMessage("Make the result layers")
startTime = float(time.time())
cx= inRaster.meanCellWidth
cy= inRaster.meanCellHeight
mx = inRaster.extent.XMin+ cx
my = inRaster.extent.YMin

llc=arcpy.Point(mx,my)
sr=inRaster.spatialReference
outRaster1 = arcpy.NumPyArrayToRaster(data_o1,llc,cx,cy,0)
arcpy.DefineProjection_management(outRaster1,sr)
arcpy.RasterToPolygon_conversion(outRaster1,outLayer,"SIMPLIFY","VALUE")
arcpy.SetParameterAsText(1,outLayer)
endTime = float(time.time())
arcpy.AddMessage(endTime-startTime)
arcpy.AddMessage("Finish")