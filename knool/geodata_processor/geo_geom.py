from osgeo import gdal,osr,ogr
import numpy as np
from PIL import Image
import zipfile
import os
import sys
import glob
import pandas as pd
import xarray as xr


def create_lineString(latlon_list): #or latlon_array
    line = ogr.Geometry(ogr.wkbLineString)
    for latlon in latlon_list:
        line.AddPoint(latlon[1], latlon[0])
    return line

def create_polygon(latlon_list): #or_latlon_array
    ring = ogr.Geometry(ogr.wkbLinearRing)
    for latlon in latlon_list:
        ring.AddPoint(latlon[1], latlon[0])
    ring.AddPoint(latlon_list[0][1], latlon_list[0][0])
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
    return poly

def create_polygons(latlon_list): #or latlon_array [n_polygon, n_point, 2(latlon)]
    multipoly = ogr.Geometry(ogr.wkbMultiPolygon)
    poly_list=[]
    for polygon in latlon_list:
        ring=ogr.Geometry(ogr.wkbLinearRing)
        [ring.AddPoint(point[1], point[0]) for point in polygon]
        ring.AddPoint(polygon[0][1], polygon[0][0])
        poly = ogr.Geometry(ogr.wkbPolygon)
        poly.AddGeometry(ring)
        poly_list.append(poly)
    return poly_list

def create_multiPolygon(latlon_list): #or latlon_array [n_polygon, n_point, 2(latlon)]
    multipoly = ogr.Geometry(ogr.wkbMultiPolygon)
    for polygon in latlon_list:
        ring=ogr.Geometry(ogr.wkbLinearRing)
        [ring.AddPoint(point[1], point[0]) for point in polygon]
        poly = ogr.Geometry(ogr.wkbPolygon)
        poly.AddGeometry(ring)
        multipoly.AddGeometry(poly)
    return multipoly

def find_intersect_polygons(poly1,poly2_list):
    # bool_list = [poly1.Intersects(poly2) for poly2 in poly2_list]
    # poly2_init_intersect=[i for (i,v) in zip(poly2_list,bool_list) if v]    
    # poly2_intersect=[poly2.Intersection(poly1) for poly2 in poly2_init_intersect]
    # poly2_init_intersect2=[polyA for polyA,polyB in zip(poly2_init_intersect,poly2_intersect) if polyB != None and polyB.Area()>0]
    # poly2_intersect2=[poly for poly in poly2_intersect if poly != None and poly.Area()>0]
    # poly2_intersect0=[]
    # for poly2 in poly2_list:
    #     print(poly2)
    #     print("hi")
    #     print(poly2.Intersection(poly1))
    #     test=poly2.Intersection(poly1)
    #     poly2_intersect0.append(test)
    poly2_intersect0=[poly2.Intersection(poly1) for poly2 in poly2_list]
    poly2_intersect_bool=[True if poly != None and poly.Area()>0 else False for poly in poly2_intersect0] 
    poly2_intersect=[poly for poly in poly2_intersect0 if poly != None and poly.Area()>0]
    return poly2_intersect, poly2_intersect_bool

