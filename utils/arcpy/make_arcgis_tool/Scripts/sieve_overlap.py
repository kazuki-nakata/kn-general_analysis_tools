#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import numpy as np
import os
import arcpy
import respack
import time
import shutil
import glob

#-----------initialize------------
fort_path=os.path.split(os.path.split(respack.__file__)[0])[0]

arcpy.gp.overwriteOutput=True

inDir = arcpy.GetParameterAsText(0)
inMaskDir = arcpy.GetParameterAsText(1)
outLayer = arcpy.GetParameterAsText(2)

# inDir = r"D:\JR\2019\ツールテスト\絞り込み\結果"
# outLayer =r"D:\RESTEC_ARCGIS_TOOLS_DEV\Scripts\tmp\test.shp"
# inMaskDir = r"D:\JR\2019\ツールテスト\絞り込み\マスク"


path=os.path.dirname(os.path.abspath(__file__))
path=path + os.sep + "tmp"

if os.path.isdir(path):
    shutil.rmtree(path)
    os.mkdir(path)
else:
    os.mkdir(path)    
    
arcpy.AddMessage(inDir)
arcpy.AddMessage(inMaskDir)
arcpy.AddMessage(outLayer)


maskShp = glob.glob(inMaskDir + os.sep + u"*.shp")
inputShp = glob.glob(inDir + os.sep + u"*.shp")

mergeShp= path + os.sep + "merge.shp"
mergeShp2= path + os.sep + "merge2.shp"
jointShp= path + os.sep + "joint.shp"

arcpy.Merge_management(inputShp, mergeShp)
arcpy.Merge_management(maskShp, mergeShp2)
arcpy.SpatialJoin_analysis(mergeShp, mergeShp2, jointShp, "#","#","#","WITHIN_A_DISTANCE",5)

with arcpy.da.UpdateCursor(jointShp, "Join_Count") as cursor:
    for row in cursor:
        if row[0] == 0:
            cursor.deleteRow()

arcpy.Copy_management(jointShp, outLayer)


