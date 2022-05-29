#!/usr/bin/env python
import sys
import numpy as np
import os
import respack
import arcpy
import scipy
import time
import ConfigParser
#-----------initialize------------
fort_path=os.path.split(os.path.split(respack.__file__)[0])[0]

tool_path=os.path.split(fort_path)[0]

arcpy.gp.overwriteOutput=True

inLayer = arcpy.GetParameterAsText(0)
inLayer2 = arcpy.GetParameterAsText(1)
outLayer = arcpy.GetParameterAsText(2)


arcpy.AddMessage("inLayer = {}".format(inLayer))
arcpy.AddMessage("inLayer = {}".format(inLayer2))
arcpy.AddMessage("outLayer = {}".format(outLayer))


