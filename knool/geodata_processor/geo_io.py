from osgeo import gdal,osr,ogr
import numpy as np
from PIL import Image
import zipfile
import os
import sys
import glob
import pandas as pd
import xarray as xr


def make_raster_from_array(data, filepath, pixel_x, pixel_y, num_band, dtype, no_data, file_type, geomode="proj", geotrans=None, geoproj=None, gcps=None, gcpsrc=None):
    
    data_mask=np.where(data == no_data, np.NaN, data)
    
    driver = gdal.GetDriverByName(file_type)    
        
    ds = driver.Create(filepath, pixel_x, pixel_y, num_band, dtype)
    
    if num_band == 1:
        band=ds.GetRasterBand(1)
        band.WriteArray(data_mask)
        band.SetNoDataValue(np.NaN)
        band.FlushCache()
    else:
        for i in range(0,num_band):
            band=ds.GetRasterBand(i+1)
            band.WriteArray(data_mask[i,:,:])
            band.SetNoDataValue(np.NaN)
            band.FlushCache()

    
    if geomode=="proj":
        ds.SetGeoTransform(geotrans)
        ds.SetProjection(geoproj)
    elif geomode == "gcp":
        ds.SetGCPs(gcps,gcpsrc)
        
    ds.FlushCache()
    return ds

def make_raster_from_array_and_prop(data, filepath, dtype, no_data, prop, file_type):
    
    if(len(data.shape) == 2):
        num_band=1
    else:
        num_band=data.shape[0]
    
    pixel_x = prop[0]
    pixel_y = prop[1]
    geotrans = prop[2]
    geoproj = prop[3]
    
    make_raster_from_array(data, filepath, pixel_x, pixel_y, num_band, dtype, no_data, file_type,"proj",geotrans=geotrans, geoproj=geoproj)

def make_raster_with_gcps_from_array(data, filepath, dtype, no_data, prop, file_type):
    
    if(len(data.shape) == 2):
        num_band=1
    else:
        num_band=data.shape[0]

    pixel_x = prop[0]
    pixel_y = prop[1]
    gcpsrc = prop[2]
    gcps = prop[3]
    
    for gcp in gcps:
        lon = gcp.GCPX
        if lon < 0:
            gcp.GCPX=lon+360.

    make_raster_from_array(data, filepath, pixel_x, pixel_y, num_band, dtype, no_data, file_type, "gcp", gcps=gcps, gcpsrc=gcpsrc)

def make_geoproj(epsg):
    geoproj = osr.SpatialReference()
    geoproj.ImportFromEPSG(epsg)
    geoproj_wkt = geoproj.ExportToWkt()
    return geoproj_wkt

def make_south_nsidc_geoinfo(res):
    geotrans = ((-3950)*1000, res*1000, 0.0, (4350)*1000, 0.0, -res*1000)
    geoproj = make_geoproj(3412)
    return geotrans, geoproj

def make_north_nsidc_geoinfo(res):
    geotrans = ((-3850)*1000, res*1000, 0.0, (5850)*1000, 0.0, -res*1000)  
    geoproj = make_geoproj(3411)
    return geotrans, geoproj

def open_generic_binary(infile,in_dtype,out_dtype,band,length,width):
    with open(infile,mode='rb') as f:
        data= np.fromfile(f, dtype=in_dtype,sep='').astype(out_dtype).reshape(band,length,width)
    if byte_order=='big':
        data= data.byteswap()
    return data

def open_hdf(infile,in_dtype,out_dtype,band,length,width):
    ds=gdal.Open(hdfpath, gdal.GA_ReadOnly)
    for i,val in enumerate(ds.GetSubDatasets()):
        name=val#[0].split(':')[4]
        ds[name]=gdal.Open(ds.GetSubDatasets()[i][0], gdal.GA_ReadOnly)
    return ds

def make_south_nsidc_raster_from_polygon(infile, outfile, width, length, band, file_type, in_dtype, byte_order, out_dtype, res,no_data):
    source_ds = ogr.Open(infile)
    source_layer = source_ds.GetLayer()
    
    geotrans,geoproj = make_south_nsidc_geoinfo(res)

    driver = gdal.GetDriverByName(file_type) 

    ds = driver.Create(outfile, width, length, band, out_dtype[1])
    ds.SetGeoTransform(geotrans)
    ds.SetProjection(geoproj)
    
    band=ds.GetRasterBand(1)
    band.SetNoDataValue(0)
    band.FlushCache()
    
    gdal.RasterizeLayer(ds, [1], source_layer,options=["ATTRIBUTE=id",'ALL_TOUCHED=TRUE']) #, burn_values=[1]
#    gdal.RasterizeLayer(ds, [1], source_layer,burn_values=[1])  
    ds.FlushCache()

def create_line_string(infile,latlon_list,epsg=4326):
    line = ogr.Geometry(ogr.wkbLineString)
    for latlon in latlon_list:
        line.AddPoint(latlon[0], latlon[1])

    ext=os.path.splitext(infile)[::-1][0]

    if ext == ".shp": 
        ftype="ESRI Shapefile"
    elif ext == ".kml":
        ftype="KML"
    else:
        print("The file type is not supported currently.")
        return
        
    driver = ogr.GetDriverByName(ftype)
    ds = driver.CreateDataSource(infile)
    srs =  osr.SpatialReference()
    srs.ImportFromEPSG(epsg)
    
    # create one layer 
    layer = ds.CreateLayer("line", srs, ogr.wkbLineString)
    # Add an ID field
    idField = ogr.FieldDefn("id", ogr.OFTInteger)
    layer.CreateField(idField)
    # Create the feature and set values
    featureDefn = layer.GetLayerDefn()
    feature = ogr.Feature(featureDefn)
    feature.SetGeometry(line)
    feature.SetField("id", 1)
    layer.CreateFeature(feature)
    feature = None
    # Save and close DataSource
    ds = None

def create_polygon(infile,latlon_list,epsg=4326):
    ring = ogr.Geometry(ogr.wkbLinearRing)
    for latlon in latlon_list:
        ring.AddPoint(latlon[0], latlon[1])
    ring.AddPoint(latlon_list[0][0], latlon_list[0][1])
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
    
    ext=os.path.splitext(infile)[::-1][0]

    if ext == ".shp": 
        ftype="ESRI Shapefile"
    elif ext == ".kml":
        ftype="KML"
    else:
        print("The file type is not supported currently.")
        return
        
    driver = ogr.GetDriverByName(ftype)
    ds = driver.CreateDataSource(infile)
    srs =  osr.SpatialReference()
    srs.ImportFromEPSG(epsg)
    
    # create one layer 
    layer = ds.CreateLayer("polygon", srs, ogr.wkbPolygon)
    # Add an ID field
    idField = ogr.FieldDefn("id", ogr.OFTInteger)
    layer.CreateField(idField)
    # Create the feature and set values
    featureDefn = layer.GetLayerDefn()
    feature = ogr.Feature(featureDefn)
    feature.SetGeometry(poly)
    feature.SetField("id", 1)
    layer.CreateFeature(feature)
    feature = None
    # Save and close DataSource
    ds = None    