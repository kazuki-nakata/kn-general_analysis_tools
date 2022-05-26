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

def create_polygon(latlon_list):
    ring = ogr.Geometry(ogr.wkbLinearRing)
    for latlon in latlon_list:
        ring.AddPoint(latlon[0], latlon[1])
    ring.AddPoint(latlon_list[0][0], latlon_list[0][1])
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
    return poly

def intersect_polygons(poly1,poly2_list):
    poly_result_list=[]
    for poly2 in poly2_list:
        if (poly1.Intersects(poly2)) # returns True
            poly_result_list.append(poly1.Intersection(poly2)) # returns None
    # intersect_geom = polyA.Intersection(polyB)
    # if(intersect_geom is not None and intersect_geom.Area()>0):
    #     return intersect_geom