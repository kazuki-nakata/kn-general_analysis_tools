from osgeo import gdal,osr,ogr
import numpy as np
from PIL import Image
import pandas as pd
import xarray as xr
from math import sin, cos, sqrt, atan2, radians, atan

def get_wgs84_wkt():
    wgs84_wkt = """
    GEOGCS["WGS 84",
        DATUM["WGS_1984",
            SPHEROID["WGS 84",6378137,298.257223563,
                AUTHORITY["EPSG","7030"]],
            AUTHORITY["EPSG","6326"]],
        PRIMEM["Greenwich",0,
            AUTHORITY["EPSG","8901"]],
        UNIT["degree",0.01745329251994328,
            AUTHORITY["EPSG","9122"]],
        AUTHORITY["EPSG","4326"]]"""
    return wgs84_wkt

def get_latlons_from_raster(raster,interval):
    # create the new coordinate system
    old_cs= osr.SpatialReference()
    old_cs.ImportFromWkt(raster.GetProjectionRef())
    wgs84_wkt = get_wgs84_wkt()
    new_cs = osr.SpatialReference()
    new_cs.ImportFromWkt(wgs84_wkt)

    transform = osr.CoordinateTransformation(old_cs,new_cs) 

    width = raster.RasterXSize
    height = raster.RasterYSize
    gt = raster.GetGeoTransform()
    minx = gt[0] + width*gt[1] + height*gt[2]
    miny = gt[3] + width*gt[4] + height*gt[5]
    xorder_vector=np.arange(0,width,interval)
    yorder_vector=np.arange(0,height,interval)
    
    loc_x=np.array([ gt[0] + xorder_vector*gt[1] + yorder*gt[2] for yorder in yorder_vector])
    loc_y=np.array([ gt[3] + xorder_vector*gt[4] + yorder*gt[5] for yorder in yorder_vector])
    cols=loc_x.shape[0]
    rows=loc_x.shape[1]
    loc_p=np.array([loc_x.reshape(cols*rows),loc_y.reshape(cols*rows)]).T #.tolist()
    latlon = transform.TransformPoints(loc_p) 
    latlon = np.array(latlon)[:,0:2].reshape(cols,rows,2).transpose((2,0,1))
    return latlon

def get_latlonval_vectors(lat_array,lon_array,val_array,mask):
    cols,rows=lat_array.shape()
    lat_vector=lat_array.shape(cols*rows)
    lon_vector=lon_array.shape(cols*rows)
    val_vector=val_array.shape(cols*rows)
    return lat_vector,lon_vector,val_vector

def get_raster_extent(raster):
    ext = []
    gt = raster.GetGeoTransform()

    xmin = gt[0]
    ymin = gt[3] + (gt[5] * raster.RasterYSize)
    xmax = gt[0] + (gt[1] * raster.RasterXSize)
    ymax = gt[3]

    ext.append(xmin)
    ext.append(ymin)
    ext.append(xmax)
    ext.append(ymax)

    return ext

def get_extent_from_corners(corners):
    ext = []
    ext.append(min([corners[0][0],corners[1][0],corners[2][0],corners[3][0]]))
    ext.append(min([corners[0][1],corners[1][1],corners[2][1],corners[3][1]]))
    ext.append(max([corners[0][0],corners[1][0],corners[2][0],corners[3][0]]))
    ext.append(max([corners[0][1],corners[1][1],corners[2][1],corners[3][1]]))
    return ext

def get_raster_corners(raster):
    ext = []
    gt = raster.GetGeoTransform()
    xarr = [0, raster.RasterXSize]
    yarr = [0, raster.RasterYSize]

    for px in xarr:
        for py in yarr:
            x = gt[0] + (px * gt[1]) + (py * gt[2])
            y = gt[3] + (px * gt[4]) + (py * gt[5])
            ext.append([x, y])
        yarr.reverse()
    return ext

def get_property_from_raster_with_proj(raster):
    prop = []
    prop.append(raster.RasterXSize)
    prop.append(raster.RasterYSize)
    prop.append(raster.GetGeoTransform())
    prop.append(raster.GetProjection())
    prop.append(osr.SpatialReference(wkt=raster.GetProjection()).GetAttrValue('AUTHORITY',1))
    prop.append(raster.GetDriver())
    return prop

def get_property_from_raster_with_gcps(raster):
    prop = []
    prop.append(raster.RasterXSize)
    prop.append(raster.RasterYSize)
    prop.append(raster.GetGCPSpatialRef())
    prop.append(raster.GetGCPs())
    prop.append(raster.GetDriver())

    if not prop[2]:
        prop[2] = """
                    GEOGCS["WGS 84",
                    DATUM["WGS_1984",
                    SPHEROID["WGS 84",6378137,298.257223563,
                        AUTHORITY["EPSG","7030"]],
                    AUTHORITY["EPSG","6326"]],
                    PRIMEM["Greenwich",0,
                        AUTHORITY["EPSG","8901"]],
                    UNIT["degree",0.01745329251994328,
                        AUTHORITY["EPSG","9122"]],
                        AUTHORITY["EPSG","4326"]]"""
    return prop

def calc_distances(lons1,lats1,lons2,lats2):
    R = 6373.0
    lon1=np.radians(lons1)
    lat1=np.radians(lats1)
    lon2=np.radians(lons2)
    lat2=np.radians(lats2)
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    distance = R * c
    return distance

def convert_lla_to_ecef(lat0, lon0, alt, a=6378137.0, b=6356752.314245):
    f = (a - b) / a
    
    lon=np.radians(lon0)
    lat=np.radians(lat0)
    
    rad = np.float64(a)        # Radius of the Earth (in meters)

    cosLat = np.cos(lat)
    sinLat = np.sin(lat)
    FF     = (1.0-f)**2
    C      = 1/np.sqrt(cosLat**2 + FF * sinLat**2)
    S      = C * FF
    x = (rad * C + alt)*cosLat * np.cos(lon)
    y = (rad * C + alt)*cosLat * np.sin(lon)
    z = (rad * S + alt)*sinLat
    return x, y, z


def convert_ecef_to_lla(x, y, z, a=6378137.0, b=6356752.314245):
    f = (a - b) / a

    e_sq = f * (2 - f)                       
    eps = e_sq / (1.0 - e_sq)
    p = np.sqrt(x * x + y * y)
    q = np.arctan2((z * a), (p * b))

    sin_q = np.sin(q)
    cos_q = np.cos(q)

    sin_q_3 = sin_q * sin_q * sin_q
    cos_q_3 = cos_q * cos_q * cos_q

    phi = np.arctan2((z + eps * b * sin_q_3), (p - e_sq * a * cos_q_3))
    lam = np.arctan2(y, x)

    v = a / np.sqrt(1.0 - e_sq * np.sin(phi) * np.sin(phi))
    h   = (p / np.cos(phi)) - v

    lat = np.degrees(phi)
    lon = np.degrees(lam)

    return lat, lon, h

def convert_ecef_to_enu(x0, y0, z0, lat0, lon0, h0, x, y, z):
    
    if type(lat0).__module__ != "numpy":
        proc=2
    elif isinstance(lat0.shape,tuple):
        proc=2
    else:
        proc=1

    if proc==1:
        px=np.radians([90 for i in range(x0.shape[0])])
        py=np.radians(90-lat0)
        pz=np.radians(lon0)
        pz2=np.radians([90 for i in range(x0.shape[0])])
        p1=np.array([x,y,z])
        p0=np.array([x0,y0,z0])
        dp=p1-p0
        
        Ry = np.array([[np.cos(py), np.zeros(x0.shape), -np.sin(py)],
                    [np.zeros(x0.shape), np.ones(x0.shape), np.zeros(x0.shape)],
                    [np.sin(py), np.zeros(x0.shape), np.cos(py)]]).transpose([2,0,1])
        
        Rz = np.array([[np.cos(pz), np.sin(pz), np.zeros(x0.shape)],
                    [-np.sin(pz), np.cos(pz), np.zeros(x0.shape)],
                    [np.zeros(x0.shape), np.zeros(x0.shape), np.ones(x0.shape)]]).transpose([2,0,1])
        
        Rz2 = np.array([[np.cos(pz2), np.sin(pz2), np.zeros(x0.shape)],
                    [-np.sin(pz2), np.cos(pz2), np.zeros(x0.shape)],
                    [np.zeros(x0.shape), np.zeros(x0.shape), np.ones(x0.shape)]]).transpose([2,0,1])
        
        R = np.matmul(np.matmul(Rz2,Ry),Rz)
        print(R.shape,dp.shape)
        result= R.dot(dp)
        x=np.diag(result[:,0,:])
        y=np.diag(result[:,1,:])
        z=np.diag(result[:,2,:])

    else:
        px=np.radians(90)
        py=np.radians(90-lat0)
        pz=np.radians(lon0)
        pz2=np.radians(90)
        p1=np.array([x,y,z])
        p0=np.array([x0,y0,z0])

        Ry = np.array([[np.cos(py), 0, -np.sin(py)],
                    [0, 1, 0],
                    [np.sin(py), 0, np.cos(py)]])
        
        Rz = np.array([[np.cos(pz), np.sin(pz), 0],
                    [-np.sin(pz), np.cos(pz), 0],
                    [0, 0, 1]])
        
        Rz2 = np.array([[np.cos(pz2), np.sin(pz2), 0],
                    [-np.sin(pz2), np.cos(pz2), 0],
                    [0, 0, 1]])

        R = Rz2.dot(Ry).dot(Rz)
        x,y,z = R.dot(p1-p0)       
        
    return x,y,z


def convert_enu_to_ecef(x0, y0, z0, lat0, lon0, h0, x, y, z):
    py=np.radians(90-lat0)
    pz=np.radians(lon0)
    pz2=np.radians(90)
    p0=np.array([x0,y0,z0])
    p2=np.array([x,y,z])

    Ry = np.array([[np.cos(py), 0, -np.sin(py)],
                   [0, 1, 0],
                   [np.sin(py), 0, np.cos(py)]])
    
    Rz = np.array([[np.cos(pz), np.sin(pz), 0],
                   [-np.sin(pz), np.cos(pz), 0],
                   [0, 0, 1]])
    
    Rz2 = np.array([[np.cos(pz2), np.sin(pz2), 0],
                   [-np.sin(pz2), np.cos(pz2), 0],
                   [0, 0, 1]])
    R = np.linalg.inv(Rz2.dot(Ry).dot(Rz))
    p1 = R.dot(p2) + p0

    return p1

def calc_line_buffer_point(lat0,lon0,h0,lat,lon,h,distance,ori="right"):
    R=6373000.0
    x0, y0, z0 = convert_lla_to_ecef(lat0, lon0, h0, a=R, b=R)
    x, y, z = convert_lla_to_ecef(lat, lon, h, a=R, b=R)    
    p, q, r = convert_ecef_to_enu(x0, y0, z0, lat0, lon0 , h0, x, y, z)
    
    sign=1
    if ori=="left" : sign=-1
    
    theta=distance/R
    tan_th=np.tan(theta)
    distance2=R*tan_th
    coef=1/(q**2/(p**2)+1)
    q2 = np.sqrt(distance2**2 * coef)
    q2= sign*q2
    p2 = -q*q2/p
    r2 = 0

    buff_ecef=convert_enu_to_ecef(x0, y0, z0, lat0, lon0 , h0, p2, q2, r2)
    lat,lon,h=convert_ecef_to_lla(*buff_ecef, a=R, b=R)
    
    return lat, lon, h


def count_shp_in_polygon(shape, polygon):
    poly = ogr.Open(polygon, 1)
    poly_lyr = poly.GetLayer(0)
    new_field = ogr.FieldDefn("Count", ogr.OFTReal)
    new_field.SetWidth(10)
    poly_lyr.CreateField(new_field)	
	
    shp = ogr.Open(shape, 0)
    shp_lyr = shp.GetLayer(0)
	
    for feature in poly_lyr:
        ext = feature.GetGeometryRef()
        shp_lyr.SetSpatialFilter(ext)
        count = shp_lyr.GetFeatureCount()
        feature.SetField("Count", count)
        poly_lyr.SetFeature(feature)
        feature.Destroy()
        shp_lyr.SetSpatialFilter(None)
		
    poly.Destroy()
    shp.Destroy()
