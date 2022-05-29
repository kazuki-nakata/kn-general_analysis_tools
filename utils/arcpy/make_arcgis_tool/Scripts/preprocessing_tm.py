#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import numpy as np
import os
import respack
import arcpy
import time
import shutil
#-----------initialize------------
fort_path=os.path.split(os.path.split(respack.__file__)[0])[0]

tool_path=os.path.split(fort_path)[0]


arcpy.gp.overwriteOutput=True


inLayer = arcpy.GetParameterAsText(0)
inMaskLayer = arcpy.GetParameterAsText(1)
outLayer = arcpy.GetParameterAsText(2)

#inLayer = r"D:\JR\2019\ORD3803_SPOT\result\IMG_SPOT6_PMS_201810290113215_ORT_SPOT6_20190814_090407zy88sq3h9mi3_1_R1C1.TIF"
#outLayer =r"D:\RESTEC_ARCGIS_TOOLS_DEV\Scripts\tmp\test.tif"
#inMaskLayer = r"D:\JR\2019\ORD3803_SPOT\result\ROI3\55_10_22_5\mask_srs.shp"
# inLayer = r"D:\JR\2019\ツールテスト\マスク\IMG_SPOT6_PMS_201810290113215_ORT_SPOT6_20190814_090834zzujlunkgpa5_1_R1C1.TIF"
# outLayer =r"D:\JR\2019\ツールテスト\マスク\result.tif"
# inMaskLayer = r"D:\JR\2019\ツールテスト\マスク\LC_N35E138.tif"
#inMaskLayer = r"D:\JR\2019\ORD3803_SPOT\LC_N35E138.tif"
ext = os.path.splitext(inMaskLayer)[1]

desc = arcpy.Describe(inLayer)
layerPath = desc.catalogPath
descPath = desc.path

inRaster = arcpy.Raster(inLayer)
inArray = arcpy.RasterToNumPyArray(inRaster)


#--------calculation------------
arcpy.AddMessage("Calculation")


#------test----------------------
cx= inRaster.meanCellWidth
cy= inRaster.meanCellHeight
mx = inRaster.extent.XMin
my = inRaster.extent.YMin
llc=arcpy.Point(mx,my)
sr=inRaster.spatialReference



path=os.path.dirname(os.path.abspath(__file__))
path=path + os.sep + "tmp"

if os.path.isdir(path):
    shutil.rmtree(path)
    os.mkdir(path)
else:
    os.mkdir(path) 

arcpy.AddMessage(path)
arcpy.CreateFeatureclass_management(path,"image_range.shp","POLYGON" )

myExtent = inRaster.extent
xmax = myExtent.XMax
xmin = myExtent.XMin
ymax = myExtent.YMax
ymin = myExtent.YMin
coordinates = [(xmin, ymin),
               (xmin, ymax),
               (xmax, ymax),
               (xmax, ymin)]
image_range = arcpy.management.CreateFeatureclass(path,
                                             "image_range.shp",
                                              "POLYGON",
                                              "",
                                              "",
                                              "",
                                              sr) # define projection
# Create feature class
outPolyExtent= image_range[0]

# Use Insert cursor to add new geometry to feature class Write feature to new feature class
with arcpy.da.InsertCursor(outPolyExtent, ['SHAPE@']) as cursor:
    cursor.insertRow([coordinates])

# # Return the Spatial Analysis extension 
arcpy.CheckInExtension("Spatial")


arcpy.CreateFileGDB_management(path, "temporary")
arcpy.env.workspace =  path + os.sep + "temporary.gdb"
#arcpy.env.compression = "NONE"
#arcpy.env.pyramid = "NONE"
arcpy.env.tileSize = "1 1"

arcpy.env.cellSize = inRaster
arcpy.env.outputCoordinateSystem = inRaster
arcpy.env.extent = inRaster
arcpy.env.snapRaster = inRaster

arcpy.AddMessage(xmin)
print ext
if ext == ".shp" or ext == "":
#---------------Shape data mask--------------------------------------
    arcpy.AddMessage("Shape data mask")
    clipFeature = path + os.sep + "clip.shp"
    projFeature = path + os.sep + "clip_proj.shp"    
    imageRange = path + os.sep + "image_range.shp"
    field = "MASK"
    maskRaster = path + os.sep + "temporary.gdb" + os.sep + "mask"
    compRaster = path + os.sep + "temporary.gdb" + os.sep + "comp_raster"

    cellSize = 1.5
    
    arcpy.Project_management(inMaskLayer, projFeature, sr)
    
    arcpy.Clip_analysis(projFeature, imageRange, clipFeature)

    arcpy.AddField_management(clipFeature, "MASK", "LONG")
    with arcpy.da.UpdateCursor(clipFeature, "MASK") as cursor:
        for row in cursor:
            row[0] = 5
            cursor.updateRow(row)

    
    arcpy.FeatureToRaster_conversion(clipFeature, field, maskRaster, cellSize)  
    arcpy.CompositeBands_management(inLayer + ";" + maskRaster,compRaster)
    arcpy.CopyRaster_management(compRaster, outLayer, "", "", "4095")

elif ext == ".tif" or ext == ".TIF" or ext == ".tiff" or ext == ".TIFF":
#---------------Shape data mask--------------------------------------  
    arcpy.AddMessage("Raster data mask")
    imageRange = path + os.sep + "image_range.shp"
    maskRaster = path + os.sep + "temporary.gdb" + os.sep + "mask"
    maskRaster = path  + os.sep + "mask.tif"    
    resampleRaster = path + os.sep + "mask_resample.tif"
    compRaster = path + os.sep + "temporary.gdb" + os.sep + "comp_raster" 
    
    arcpy.Clip_management(inMaskLayer,"", maskRaster, imageRange, "0", "ClippingGeometry")    
#    cellSize = (cx,cy)
    arcpy.ProjectRaster_management(maskRaster,resampleRaster,sr,"NEAREST",1.5)
#    arcpy.Resample_management(maskRaster, resampleRaster, cellSize, "NEAREST")
    arcpy.CompositeBands_management(inLayer + ";" + resampleRaster,compRaster)
#    arcpy.CompositeBands_management(resampleRaster + ";" + inLayer,compRaster)
    
    arcpy.CopyRaster_management(compRaster, outLayer, "", "", "")
    arcpy.AddMessage(xmin)
else:
    arcpy.AddMessage("Select a ESRI Shape file (.shp) or TIFF image (.tif, .TIF, .tiff, .TIFF) as mask data")
