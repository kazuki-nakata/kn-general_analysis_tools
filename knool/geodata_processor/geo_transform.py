from osgeo import gdal, osr, ogr
import os
from . import ogr2ogr
from . import gdal_merge as merge
from . import geo_info

OGRTypes = {int: ogr.OFTInteger, str: ogr.OFTString, float: ogr.OFTReal}


def reverse_clip(in_maskFile, out_maskFile):
    #    options = ["", "-f", "ESRI Shapefile", "-sql",
    #    ''"SELECT ST_Difference(a.geometry, b.geometry) AS geometry FROM a,'b.shp'.b"'',
    #    wkt_pr, outmaskFile, maskFile]
    #    ogr2ogr difference.shp a.shp -dialect SQLite -sql
    geoproj = osr.SpatialReference()
    geoproj.ImportFromEPSG(3411)
    wkt_pr = geoproj.ExportToWkt()
    options = ["", "-f", "ESRI Shapefile", "-t_srs", wkt_pr, out_maskFile, in_maskFile]
    ogr2ogr.main(options)


def reproject(ds, outfile="/vsimem/output.tif", epsg_str="EPSG:4326", NODATA_VALUE=0):
    output_ds = gdal.Warp(outfile, ds, dstSRS=epsg_str, resampleAlg="bilinear", dstNodata=NODATA_VALUE)
    return output_ds


def adjust_shape(ds, ds2, outfile="/vsimem/output.tif"):
    cols = ds2.RasterXsize
    rows = ds2.RasterYSize
    output_ds = gdal.Warp(outfile, ds, width=cols, height=rows)
    return output_ds


def clip_rectangle(ds, minX, minY, maxX, maxY, outfile="/vsimem/output.tif"):
    output_ds = gdal.Warp(outfile, ds, outputBounds=(minX, minY, maxX, maxY))
    return output_ds


def warp(ds, outfile, epsg_str, xres, yres, minX, minY, maxX, maxY, resample):
    output_ds = gdal.Warp(
        outfile, ds, dstSRS=epsg_str, xRes=xres, yRes=yres, outputBounds=(minX, minY, maxX, maxY), resampleAlg=resample
    )
    return output_ds


# def edit_tifftag(ds, tag_name, tag_value, outfile="/vsimem/output.tif"):
#     # For th edetails, see gdal official HP https://gdal.org/drivers/raster/gtiff.html
#     # TFW=YES :
#     # RPB=YES : RPC 情報が利用可能な場合
#     # TILED=YES : デフォルトでは、ストライプ化された TIFF ファイルが作成されます。このオプションは、タイル化された TIFF ファイルの作成を強制するために使用できます。
#     # BLOCKXSIZE=n : タイル幅を設定します。デフォルトは 256 です。
#     # BLOCKYSIZE=n : タイルまたはストリップの高さを設定します。タイルの高さのデフォルトは 256 で、ストリップの高さのデフォルトは、1 つのストリップが 8K 以下になるような値です。
#     # COMPRESS=[JPEG/LZW/PACKBITS/DEFLATE/CCITTRLE/CCITTFAX3/CCITTFAX4/LZMA/ZSTD/LERC/LERC_DEFLATE/LERC_ZSTD/WEBP/JXL/NONE]
#     ds = gdal.Translate(
#         outfile,
#         ds,
#     )
#     return ds


def collocate(outfile, infile1, infile2, ulx, uly, lrx, lry):
    # fmt: off
    options = ["", "-o", outfile, "-separate", "-n", "0", "-ul_lr",
               str(ulx), str(uly), str(lrx), str(lry), infile1, infile2]
    # fmt: on
    print(options)
    merge.main(options)


def mosaic(infile_list, outfile):
    options = ["", "-o", outfile]
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
    gdal.Polygonize(srcband, srcband, dst_layer, -1, ["8CONNECTED=8"])


def get_vector_from_raster(raster):
    srs = osr.SpatialReference()
    srs.ImportFromWkt(raster.GetProjection())
    srcband = raster.GetRasterBand(1)
    drv = ogr.GetDriverByName("GPKG")
    dst_ds = drv.CreateDataSource(r"/vsimem/output.gpkg")
    dst_layer = dst_ds.CreateLayer("output", srs=srs)
    gdal.Polygonize(srcband, srcband, dst_layer, -1, ["8CONNECTED=8"])
    vector = ogr.Open("/vsimem/output.gpkg")
    return vector


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


def rasterize_by_raster_with_gcps(source_ds, raster_ds, outfile, file_type="GTiff", out_dtype=gdal.GDT_Float32):

    source_layer = source_ds.GetLayer()

    prop = geo_info.get_property_from_raster_with_gcps(raster_ds)

    num_band = 1
    pixel_x = prop[0]
    pixel_y = prop[1]
    gcpsrc = prop[2]
    gcps = prop[3]

    for gcp in gcps:
        lon = gcp.GCPX
        if lon < 0:
            gcp.GCPX = lon + 360.0

    driver = gdal.GetDriverByName(file_type)
    ds = driver.Create(outfile, pixel_x, pixel_y, num_band, out_dtype)
    band = ds.GetRasterBand(1)
    band.SetNoDataValue(0)
    band.FlushCache()

    #    ds.SetProjection(geoproj)
    ds.SetGCPs(gcps, gcpsrc)

    gdal.RasterizeLayer(ds, [1], source_layer, options=["ATTRIBUTE=id", "ALL_TOUCHED=TRUE"])  # , burn_values=[1]
    #    gdal.RasterizeLayer(ds, [1], source_layer,burn_values=[1])
    ds.FlushCache()


def reproject_shp(source_shp, target_shp, t_srs):
    s_ds = ogr.Open(source_shp, 0)
    s_lyr = s_ds.GetLayer(0)
    s_srs = s_lyr.GetSpatialRef()

    coord_trans = osr.CoordinateTransformation(s_srs, t_srs)

    t_driver = ogr.GetDriverByName("ESRI Shapefile")
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
    file = open(f_name + ".prj", "w")
    file.write(t_srs.ExportToWkt())
    file.close()

    s_ds.Destroy()
    t_ds.Destroy()
