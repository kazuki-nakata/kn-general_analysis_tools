from osgeo import gdal,osr,ogr
import numpy as np
import os
import ogr2ogr
import gdal_merge as merge


def reverse_clip(in_maskFile,out_maskFile):
#    options = ["", "-f", "ESRI Shapefile", "-sql", ''"SELECT ST_Difference(a.geometry, b.geometry) AS geometry FROM a, 'b.shp'.b"'', wkt_pr, outmaskFile, maskFile]
#    ogr2ogr difference.shp a.shp -dialect SQLite -sql 
    geoproj = osr.SpatialReference()
    geoproj.ImportFromEPSG(3411)
    wkt_pr = geoproj.ExportToWkt()
    options = ["", "-f", "ESRI Shapefile", "-t_srs", wkt_pr, out_maskFile, in_maskFile]
    ogr2ogr.main(options)

def reproject(ds,outfile='/vsimem/output.tif',epsg_str="EPSG:4326"):
    output_ds = gdal.Warp(outfile, ds, dstSRS=epsg_str,resampleAlg="bilinear")
    return output_ds

def adjust_shape(ds,ds2,outfile='/vsimem/output.tif'):
    cols = ds2.RasterXsize
    rows = ds2.RasterYSize
    output_ds = gdal.Warp(outfile, ds, width=cols, height=rows)
    return output_ds

def warp(ds,outfile, epsg_str,xres,yres,resample):
    output_ds = gdal.Warp(outfile, ds, dstSRS=epsg_str, xRes=xres, yRes=yres, resampleAlg=resample)
    return output_ds

def collocate(outfile,infile1,infile2,ulx,uly,lrx,lry):
    options = ['', '-o', outfile, '-separate', '-n', '0', '-ul_lr',
               str(ulx), str(uly), str(lrx), str(lry), infile1, infile2]
    print(options)
    merge.main(options)

def mosaic(infile_list,outfile):
    options = ['', '-o', outfile]
    for infile in infile_list:
        options.append(infile)
    merge.main(options)

def raster_to_polygon(raster, dst_layername):
    srs = osr.SpatialReference()
    srs.ImportFromWkt(raster.GetProjection())
    srcband = raster.GetRasterBand(1)
    drv = ogr.GetDriverByName("ESRI Shapefile")
    dst_ds = drv.CreateDataSource(dst_layername + ".shp")
    dst_layer = dst_ds.CreateLayer(dst_layername, srs=srs)
    gdal.Polygonize(srcband, srcband, dst_layer, 0, ["8CONNECTED=8"])

def raster_to_polygon_with_field(raster, dst_layername, field_name):
    srs = osr.SpatialReference()
    srs.ImportFromWkt(raster.GetProjection())
    srcband = raster.GetRasterBand(1)
    drv = ogr.GetDriverByName("ESRI Shapefile")
    dst_ds = drv.CreateDataSource(dst_layername + ".shp")
    dst_layer = dst_ds.CreateLayer(dst_layername, srs=srs)
    fd = ogr.FieldDefn(field_name, ogr.OFTInteger)
    dst_layer.CreateField(fd)
    dst_field = dst_layer.GetLayerDefn().GetFieldIndex(field_name)
    gdal.Polygonize(srcband, srcband, dst_layer, dst_field, ["8CONNECTED=8"])

def reproject_shp(source_shp, target_shp, t_srs):
    s_ds = ogr.Open(source_shp, 0)
    s_lyr = s_ds.GetLayer(0)
    s_srs = s_lyr.GetSpatialRef()

    coord_trans = osr.CoordinateTransformation(s_srs, t_srs)

    t_driver = ogr.GetDriverByName('ESRI Shapefile')
    t_ds = t_driver.CreateDataSource(target_shp)
    t_lyr_name = os.path.splitext(os.path.split(target_shp)[1])[0]
    t_lyr = t_ds.CreateLayer(t_lyr_name, geom_type=ogr.wkbMultiPolygon)

    s_lyr_defn = s_lyr.GetLayerDefn()
    for i in range(0, s_lyr_defn.GetFieldCount()):
        field_defn = s_lyr_defn.GetFieldDefn(i)
        t_lyr.CreateField(field_defn)

    t_lyr_defn = t_lyr.GetLayerDefn()

    for feature in s_lyr:
        geom = feature.GetGeometryRef()
        geom.Transform(coord_trans)
        t_feature = ogr.Feature(t_lyr_defn)
        t_feature.SetGeometry(geom)

        for i in range(0, t_lyr_defn.GetFieldCount()):
            t_feature.SetField(t_lyr_defn.GetFieldDefn(i).GetNameRef(), feature.GetField(i))

        t_lyr.CreateFeature(t_feature)

        t_feature.Destroy()
        feature.Destroy()

    t_srs.MorphToESRI()
    f_name = os.path.splitext(target_shp)[0]
    file = open(f_name + ".prj", 'w')
    file.write(t_srs.ExportToWkt())
    file.close()

    s_ds.Destroy()
    t_ds.Destroy()