from osgeo import gdal,osr,ogr
import numpy as np
from PIL import Image
import zipfile
import os
import sys
import glob
import pandas as pd
import xarray as xr

def create_line_string(latlon_list):
    line = ogr.Geometry(ogr.wkbLineString)
    for latlon in latlon_list:
        line.AddPoint(latlon[0], latlon[1])
    return line

def create_polygon(latlon_lis):
    ring = ogr.Geometry(ogr.wkbLinearRing)
    for latlon in latlon_list:
        ring.AddPoint(latlon[0], latlon[1])
    ring.AddPoint(latlon_list[0][0], latlon_list[0][1])
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
    return poly