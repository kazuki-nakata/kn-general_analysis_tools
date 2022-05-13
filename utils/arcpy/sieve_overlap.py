
import numpy as np
import os
import arcpy

import shutil
import glob

#-----------initialize------------
arcpy.gp.overwriteOutput=True

inDir = arcpy.GetParameterAsText(0)
inMaskDir = arcpy.GetParameterAsText(1)
outLayer = arcpy.GetParameterAsText(2)

inDir = r"indir where polygons are"
outLayer =r"outputpath"
inMaskDir = r"polygonpath for mask"


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
